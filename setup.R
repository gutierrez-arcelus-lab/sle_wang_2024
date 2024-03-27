library(tidyverse)
library(furrr)
library(glue)
library(httr)
library(jsonlite)
library(biomaRt)

select <- dplyr::select

# Create results directory
res_dir <- "./results" 
if (!file.exists(res_dir)) dir.create(res_dir)

qtl_data_dir <- "./data/eQTLcatalogue/sumstats" 
if (!file.exists(qtl_data_dir)) dir.create(ge_data_dir, recursive = TRUE)

if (!file.exists("data/gwas")) dir.create("data/gwas")

# eQTL catalogue datasets
max_pulled_rows <- 1000
URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")

r <- GET(URL, accept_json())
status_code(r)

cont <- content(r, "text", encoding = "UTF-8")

datasets <- 
    fromJSON(cont) |>
    as_tibble()

write_tsv(datasets, "./data/eQTLcatalogue/metadata.tsv")

cat_dir <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue"

qtl_files <-
    datasets |>
    filter(quant_method %in% c("ge", "exon")) |>
    mutate(f = glue("{cat_dir}/{quant_method}_all/{dataset_id}.all.tsv.gz")) |>
    mutate(f = ifelse(study_label == "Lepik_2017" & quant_method == "exon", 
		      glue("{cat_dir}/Lepik_2017/{dataset_id}.all.tsv.gz"),
		      f)) |>
    select(dataset_id, f) |>
    deframe()

dataset_ids <- names(qtl_files)

# Gene annotations
gene_annots <-
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "##", col_names = FALSE) |>
    filter(X3 == "gene")

# Extract IKZF1 summary stats from eQTL Catalogue
ikzf1_tss <- 
    gene_annots |>
    filter(X1 == "chr7") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    filter(gene_name == "IKZF1") |>
    pull(X4)

ikzf1_window <- 
    paste(ikzf1_tss + c(-1e6, 1e6), collapse = "-") |>
    {function(x) paste0("7:", x)}()

write_lines(ikzf1_window, "./data/ikzf1_coords_grch38.txt")

cmd <- 
    glue("tabix {qtl_files} {ikzf1_window} | bgzip > {qtl_data_dir}/{dataset_ids}.ikzf1.tsv.gz")

plan(multisession, workers = availableCores())
future_walk(cmd, system)
plan(sequential)

# GWAS summary statistics
## Langefeld
langefeld_chr7 <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr7.out" |>
    read_delim(comment = "#", delim = " ") |>
    filter(all_maf >= 0.002) |>
    select(rsid, 
	   pos = position, 
	   alleleA, 
	   alleleB,
	   beta = `frequentist_add_beta_1:add/sle=1`, 
	   se = frequentist_add_se_1, 
	   p = frequentist_add_wald_pvalue_1) |>
    drop_na(beta, se) 

ikzf1_risk_var <- langefeld_chr7 |> 
    filter(rsid == "rs4917014") |>
    pull(pos)

langefeld_ikzf1 <- langefeld_chr7 |>
    filter(between(pos, ikzf1_risk_var - 5e5, ikzf1_risk_var + 5e5)) |>
    add_column(chrom = "7", .before = 1) |>
    drop_na(beta, se, p) |>
    mutate(rsid = str_extract(rsid, "rs\\d+")) |>
    select(rsid, pos, alleleA, alleleB, beta, se, p)

# Update and get missing RsIDs by matching position and alleles to dbSNP
dbsnp_meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt" |>
    data.table::fread(skip = "# Sequence-Name") 

refseq_chr_id <- dbsnp_meta |>
    filter(`UCSC-style-name` == "chr7") |>
    pull(`RefSeq-Accn`)

dbsnp_ikzf1_window <- 
    paste0(refseq_chr_id, ":", paste(ikzf1_tss + c(-1e6, 1e6), collapse = "-"))

dbsnp_vcf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25.gz"

glue("tabix {dbsnp_vcf} {dbsnp_ikzf1_window}") |>
    paste("| awk '{ print $2,$3,$4,$5 }' > data/dbsnp_ikzf1.vcf") |>
    system()

dbsnp_ikzf1 <- "./data/dbsnp_ikzf1.vcf" |>  
    read_delim(delim = " ", col_names = c("pos", "rsid", "ref", "alt")) |>
    separate_rows(alt, sep = ",")

langefeld_ikzf1_dbsnp <- 
    langefeld_ikzf1 |>
    left_join(dbsnp_ikzf1, join_by(pos, alleleA == ref, alleleB == alt)) |>
    mutate(rsid = ifelse(!is.na(rsid.y), rsid.y, rsid.x)) |>
    select(rsid, pos, alleleA, alleleB, beta, se, p) |>
    filter(!is.na(rsid))

## Bentham ("GRCh38")
bentham_ikzf1 <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    read_tsv() |>
    drop_na(beta) |>
    filter(hm_chrom == 7, 
	   between(hm_pos, ikzf1_risk_var - 5e5, ikzf1_risk_var + 5e5)) |>
    select(chr = hm_chrom, rsid = hm_rsid, pos = hm_pos,
	   allele_a = hm_other_allele, allele_b = hm_effect_allele, 
	   beta, se = standard_error)

# Save GWAS data
write_tsv(langefeld_ikzf1_dbsnp, "./data/gwas/langefeld_ikzf1.tsv")
write_tsv(bentham_ikzf1, "./data/gwas/bentham_ikzf1.tsv")

# Save coords in hg19 for LD estimation
paste0("7:", ikzf1_risk_var - 1e6, "-", ikzf1_risk_var + 1e6) |>
    write_lines("./data/ikzf1_coords_grch37.txt")

# liftOver Langefeld data
# LiftOver
bed19_file <- "./data/gwas/langefeld_hg19.bed"
bed38_file <- "./data/gwas/langefeld_hg38.bed"
fail_file <- "./data/gwas/langefeld.failTolift.txt"
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"

bed19 <- 
    "./data/gwas/langefeld_ikzf1.tsv" |>
    read_tsv() |>
    mutate(chr = "chr7",
	   end = pos,
	   start = pos - 1L) |>
    select(chr, start, end, rsid)

write_tsv(bed19, bed19_file, col_names = FALSE)

command <- sprintf("liftOver %s %s %s %s", bed19_file, chain_file, bed38_file, fail_file)
system(command)


# Save KGP IDs for LD reference panel
dat <- 
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##")

europeans <- c("CEU", "TSI", "GBR", "IBS", "FIN")

dat |>
    distinct(SAMPLE_NAME, POPULATION) |>
    filter(POPULATION %in% europeans) |>
    arrange(SAMPLE_NAME) |>
    pull(SAMPLE_NAME) |>
    write_lines("./data/eur_kgp.txt")

