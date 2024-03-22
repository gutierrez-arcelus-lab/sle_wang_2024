library(tidyverse)
library(furrr)
library(coloc)

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)
bentham_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/opengwas/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##")

snp_window <- bentham_stats |> 
    filter(ID == "rs2297550") |>
    mutate(start = POS - 2.5e5L, stop = POS + 2.5e5L) |>
    select(`#CHROM`, start, stop)

bentham_stats_window <- bentham_stats |>
    inner_join(snp_window, join_by(`#CHROM`, between(POS, start, stop))) |>
    separate_rows(FORMAT:`EBI-a-GCST003156`, sep = ":") |>
    select(chr = `#CHROM`, pos = POS, rsid = ID, REF, ALT, FORMAT, stats = `EBI-a-GCST003156`) |>
    pivot_wider(names_from = FORMAT, values_from = stats) |>
    select(chr, pos, rsid = ID, ref = REF, alt = ALT, beta = ES, se = SE, logp = LP) |>
    mutate_at(vars(beta, se), as.numeric)

eqtl_cat <- "./data/eqtl_catalogue_associations.tsv.gz" |>
    read_tsv(col_types = c(chromosome = "c")) |>
    inner_join(snp_window, join_by(chromosome == `#CHROM`, between(position, start, stop))) |>
    select(dataset_id, gene_id, molecular_trait_id, chr = chromosome, pos = position, 
	   rsid, ref, alt, an, beta, maf, se, pvalue)

min_df <- 
    inner_join(bentham_stats_window, eqtl_cat,
	       join_by(chr, rsid, ref, alt),
	       suffix = c("_gwas", "_eqtl")) |>
    select(dataset_id, gene_id, molecular_trait_id, rsid, beta_gwas, se_gwas,
	   beta_eqtl, se_eqtl, an_eqtl = an, maf_eqtl = maf)


run_coloc <- function(min_df) {

    eqtl_dataset <- list(beta = min_df$beta_eqtl,
			 varbeta = min_df$se_eqtl^2,
			 N = min_df$an_eqtl/2,
			 MAF = min_df$maf_eqtl,
			 type = "quant",
			 snp = min_df$rsid)

    gwas_dataset <- list(beta = min_df$beta_gwas,
			 varbeta = min_df$se_gwas^2,
			 type = "cc",
			 snp = min_df$rsid)

    coloc_res <- coloc.abf(dataset1 = eqtl_dataset, 
			   dataset2 = gwas_dataset)

    coloc_res
}

coloc_res <- min_df |>
    unite("subset_id", c(dataset_id, gene_id, molecular_trait_id), sep = ",") |> 
    {function(x) split(x, x$subset_id)}() |>
    map(run_coloc)

coloc_results_summary <- map(coloc_res, "summary") |>
    map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
    pivot_wider(names_from = stat, values_from = value) |>
    select(id, nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf) |>
    separate(id, c("dataset_id", "gene_id", "molecular_trait_id"), sep = ",")

# annotations
annot <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v39.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "##", col_names = FALSE)

gene_info <- annot |>
    filter(X3 == "gene", X1 == "chr1") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(gene_id, gene_name)

# eQTL catalogue metadata
eqtl_meta <- read_tsv("./data/eqtl_datasets.tsv") |>
    select(-study_id, -tissue_id)

out <- left_join(coloc_results_summary, gene_info, join_by(gene_id)) |>
    left_join(eqtl_meta, join_by(dataset_id)) |>
    select(dataset_id, study = study_label, sample_group, tissue = tissue_label,
	   condition = condition_label, level = quant_method,
	   gene_id, gene_name, molecular_trait_id, everything())

write_tsv(out, "./data/coloc_results.tsv")

