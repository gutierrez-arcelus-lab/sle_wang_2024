library(tidyverse)
library(furrr)
library(coloc)
library(glue)

run_coloc <- function(min_df) {

    eqtl_dataset <- list(beta = min_df$beta_eqtl,
			 varbeta = min_df$se_eqtl^2,
			 N = min_df$an_eqtl/2,
			 MAF = min_df$maf_eqtl,
			 type = "quant",
			 snp = min_df$var_id)

    gwas_dataset <- list(beta = min_df$beta_gwas,
			 varbeta = min_df$se_gwas^2,
			 type = "cc",
			 snp = min_df$var_id)

    coloc_res <- coloc.abf(dataset1 = eqtl_dataset, 
			   dataset2 = gwas_dataset)

    coloc_res
}

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

# GWAS data
langefeld_stats <- 
    "./data/gwas/langefeld_ikzf1.tsv" |>
    read_tsv()

bentham_stats <- 
    "./data/gwas/bentham_ikzf1.tsv" |>
    read_tsv()

# eQTL Catalogue data
qtl_datasets <- 
    read_tsv("./data/eQTLcatalogue/metadata.tsv") |>
    select(study_id, dataset_id, study_label, sample_group, tissue_label, condition_label, quant_method)

qtl_files <- 
    "./data/eQTLcatalogue/sumstats" |>
    list.files(pattern = "QTD\\d+\\.ikzf1\\.tsv\\.gz", full.names = TRUE) |>
    {function(x) setNames(x, str_extract(x, "QTD\\d+"))}()

qtl_names <-
    "zcat /lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/ge_all/QTD000001.all.tsv.gz | head -n1" |>
    system(intern = TRUE) |>
    strsplit("\t") |>
    unlist() |>
    str_remove("\\r$")

plan(multisession, workers = availableCores())

qtl_stats <-
    qtl_files |>
    future_map_dfr(~read_tsv(., col_names = qtl_names) |> 
		   select(chr = chromosome, 
			  molecular_trait_id, 
			  gene_id, 
			  pos = position, 
			  rsid, ref, alt, an, beta, maf, se, pvalue) |>
		   mutate(p_fdr = p.adjust(pvalue, method = "fdr")) |>
		   group_by(molecular_trait_id) |>
		   filter(any(p_fdr <= 0.05)) |>
		   ungroup(),
		   .id = "dataset")

plan(sequential)

# Run coloc
## Langefeld
langefeld_min_df <-
    inner_join(qtl_stats, langefeld_stats,
	       join_by(rsid, ref == alleleA, alt == alleleB),
	       suffix = c("_eqtl", "_gwas")) |>
    unite("var_id", c("rsid", "ref", "alt"), sep = "-") |>
    select(dataset, gene_id, molecular_trait_id, var_id, 
	   beta_gwas, se_gwas, beta_eqtl, se_eqtl, an_eqtl = an, maf_eqtl = maf)

coloc_res_langefeld <- langefeld_min_df |>
    unite("subset_id", c(dataset, gene_id, molecular_trait_id), sep = ",") |> 
    {function(x) split(x, x$subset_id)}() |>
    map(run_coloc)

coloc_summary_langefeld <- map(coloc_res_langefeld, "summary") |>
    map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
    pivot_wider(names_from = stat, values_from = value) |>
    select(id, nsnps, pp0 = PP.H0.abf, pp1 = PP.H1.abf, pp2 = PP.H2.abf, pp3 = PP.H3.abf, pp4 = PP.H4.abf) |>
    separate(id, c("dataset_id", "gene_id", "molecular_trait_id"), sep = ",")

## Bentham
bentham_min_df <-
    inner_join(qtl_stats, bentham_stats,
	       join_by(rsid, ref == allele_a, alt == allele_b),
	       suffix = c("_eqtl", "_gwas")) |>
    unite("var_id", c("rsid", "ref", "alt"), sep = "-") |>
    select(dataset, gene_id, molecular_trait_id, var_id, 
	   beta_gwas, se_gwas, beta_eqtl, se_eqtl, an_eqtl = an, maf_eqtl = maf)

coloc_res_bentham <- bentham_min_df |>
    filter(beta_gwas != 0) |>
    unite("subset_id", c(dataset, gene_id, molecular_trait_id), sep = ",") |> 
    {function(x) split(x, x$subset_id)}() |>
    map(run_coloc)

coloc_summary_bentham <- map(coloc_res_bentham, "summary") |>
    map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
    pivot_wider(names_from = stat, values_from = value) |>
    select(id, nsnps, pp0 = PP.H0.abf, pp1 = PP.H1.abf, pp2 = PP.H2.abf, pp3 = PP.H3.abf, pp4 = PP.H4.abf) |>
    separate(id, c("dataset_id", "gene_id", "molecular_trait_id"), sep = ",")

coloc_out <-
    bind_rows("Langefeld" = coloc_summary_langefeld, 
	      "Bentham" = coloc_summary_bentham,
	      .id = "gwas_id")

# annotations
annot <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v39.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "##", col_names = FALSE)

gene_info <- 
    annot |>
    filter(X3 == "gene", X1 == "chr7") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(gene_id, gene_name)

# Save results
out <- 
    left_join(coloc_out, gene_info, join_by(gene_id)) |>
    left_join(qtl_datasets, join_by(dataset_id)) |>
    select(gwas_id, dataset_id, study_label, sample_group, 
	   tissue_label, condition_label, quant_method,
	   gene_id, gene_name, molecular_trait_id, 
	   nsnps, pp0:pp4)

write_tsv(out, "./results/coloc_results_ikzf1.tsv")

# Test
signif_colocs <- out |>
    filter(pp4 > .5) |>
    distinct(dataset_id, study_label, tissue_label, condition_label, quant_method, gene_id, gene_name, molecular_trait_id) |>
    filter(tissue_label %in% c("LCL", "Treg memory", "CD16+ monocyte"))

permuted <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/exon_permuted" |>
    list.files(full.names = TRUE) |>
    {function(x) setNames(x, str_extract(x, "QTD\\d+"))}() |>
    map_dfr(read_tsv, .id = "dataset_id") |>
    filter(chromosome == 7, between(position, 49e6, 51e6))

inner_join(permuted, select(signif_colocs, dataset_id, gene_id, gene_name),
	   join_by(dataset_id, molecular_trait_object_id == gene_id)) |>
    group_by(dataset_id) |>    
    mutate(p_fdr = p.adjust(p_beta, method = "fdr")) |>
    ungroup() |>
    print(width = Inf)

signif_colocs |> print(width = Inf)

"/lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/exon_permuted/QTD000106.permuted.tsv.gz" |>
    read_tsv() |>
    mutate(p_fdr = p.adjust(p_beta, method = "fdr")) |>
    filter(grepl("ENSG00000185811", molecular_trait_id)) |>
    print(width = Inf)

