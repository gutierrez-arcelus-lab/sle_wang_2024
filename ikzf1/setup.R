library(tidyverse)

if (!file.exists("data")) dir.create("data")

# Summary stats
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
	   p_lrt = frequentist_add_lrt_pvalue,
	   p_wald = frequentist_add_wald_pvalue_1) |>
    drop_na(beta, se) 

risk_var <- langefeld_chr7 |> 
    filter(rsid == "rs4917014") |>
    pull(pos)


#### Japanese
#japan <- 
#    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Japanese_SLE_20210508/METAL_Auto_X_with_GC.meta" |>
#    read_delim(delim = " ") |>
#    filter(grepl("^7:", MARKER))
#
#japan |> 
#extract(MARKER, c("chr", "pos"), "([^:]+):(\\d+)") |>
#mutate(pos = as.integer(pos)) |>
#filter(pos == risk_var)
#
####

langefeld_region <- langefeld_chr7 |>
    filter(between(pos, risk_var - 5e5, risk_var + 5e5)) |>
    add_column(chrom = "7", .before = 1) |>
    drop_na(beta, se, p_lrt)

# Variants with RsID
langefeld_region_rsid <- 
    langefeld_region |>
    filter(grepl("rs\\d+", rsid)) |>
    mutate(rsid = str_extract(rsid, "rs\\d+")) |>
    unite("var_id", c(chrom, pos, alleleA, alleleB), sep = "_", remove = FALSE) |>
    select(var_id, rsid, pos, alleleA, alleleB, beta, se, p_lrt, p_wald)

# Variants with no RsID
# Exclude struc variants
langefeld_region_norsid <- 
    langefeld_region |>
    filter(!grepl("rs\\d+", rsid), !grepl("CN|INS", rsid)) |>
    mutate(var_id = str_replace_all(rsid, ":", "_"),
	   rsid = NA) |>
    select(var_id, rsid, pos, alleleA, alleleB, beta, se, p_lrt, p_wald)
	   
langefeld_clean <- bind_rows(langefeld_region_rsid, langefeld_region_norsid)

# LiftOver
bed19_file <- "./data/langefeld_hg19.bed"
bed38_file <- "./data/langefeld_hg38.bed"
fail_file <- "./data/langefeld.failTolift.txt"
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"

bed19 <- 
    langefeld_clean |>
    mutate(chr = "chr7",
	   end = pos,
	   start = pos - 1L) |>
    select(chr, start, end, var_id)

write_tsv(bed19, bed19_file, col_names = FALSE)

command <- sprintf("liftOver %s %s %s %s", bed19_file, chain_file, bed38_file, fail_file)
system(command)

bed38 <- 
    read_tsv(bed38_file, col_names = c("chr", "start", "end", "var_id")) |>
    select(var_id, pos_hg38 = end)

langefeld_hg38 <- 
    left_join(langefeld_clean, bed38, join_by(var_id)) |>
    select(var_id, rsid, pos, pos_hg38, everything()) |>
    mutate(chr = "7") |>
    unite("var_id", c(chr, pos, alleleA, alleleB), remove = FALSE) |>
    select(rsid, pos = pos_hg38, alleleA, alleleB, beta:p_wald)

write_tsv(langefeld_hg38, "./data/langefeld_clean.tsv")

region <- sprintf("chr7:%s", paste(range(langefeld_hg38$pos), collapse = "-"))
write_lines(region, "./data/langefeld_region.txt") 



## dbSNP
#dbsnp_hg38 <- 
#    read_tsv("./data/dbsnp_ikzf1_hg38.vcf", comment = "##") |>
#    select(chrom = 1, pos = POS, rsid = ID, REF, ALT) |>
#    separate_rows(ALT, sep = ",")
#
# Match by rsid, lifted position, and alleles
#langefeld_dbsnp_match <- 
#    langefeld_hg38 |>
#    inner_join(dbsnp_hg38, join_by(rsid, pos, alleleA == REF, alleleB == ALT)) |>
#    select(rsid, pos, ref = alleleA, alt = alleleB, beta:p_wald)
#
# Usually, when variants match by position and alleles, 
# only a few RsIDs will be a mismatch due to changing RsIDs, RsIDs withdrawn, 
# or even RsIDs in Immunochiop that don't exist or never existed in dbSNP.
# I tried a few methods to get the correct RsIDs,
# eg, querying dbSNP website with a custom function or with rsnps::ncbi_snp_query.
# Surprisingly, matching by position and alleles, and getting the rsid in the dbSNP VCF 
# seems to be the best method.
# There are a few variants that match by lifted position and RsID, but not alleles.
# Sometimes they match by inversing alleleA and alleleB.
# Most of these seem to be a good match not considering alleles, meaning that
# alleles encoded as alleleA and alleleB correspond to ALT and REF, instead of REF and ALT.

# In the end, I chose to match by lifted position and by alleles, taking into account possibility 
# that alleles are flipped.
#langefeld_dbsnp_match <- 
#    langefeld_hg38 |>
#    inner_join(dbsnp_hg38, join_by(pos, alleleA == REF, alleleB == ALT)) |>
#    select(rsid = rsid.y, pos, ref = alleleA, alt = alleleB, beta:p_wald)
#
#langefeld_dbsnp_match_flip <- 
#    langefeld_hg38 |>
#    inner_join(dbsnp_hg38, join_by(pos, alleleA == ALT, alleleB == REF)) |>
#    select(rsid = rsid.y, pos, ref = alleleA, alt = alleleB, beta:p_wald) |>
#    mutate(beta = beta * -1)
#
#langefeld_clean_final <- bind_rows(langefeld_dbsnp_match, langefeld_dbsnp_match_flip)  
#
#write_tsv(langefeld_clean_final, "./data/langefeld_clean.tsv")
#
#region <- sprintf("chr7:%s", paste(range(langefeld_clean_final$pos), collapse = "-"))
#write_lines(region, "./data/langefeld_region.txt") 
#
#
# 1000 Genomes data to serve as a reference panel
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

chinese <- c("CDX", "CHB", "CHS")

dat |>
    distinct(SAMPLE_NAME, POPULATION) |>
    filter(POPULATION %in% chinese) |>
    arrange(SAMPLE_NAME) |>
    pull(SAMPLE_NAME) |>
    write_lines("./data/chi_kgp.txt")

