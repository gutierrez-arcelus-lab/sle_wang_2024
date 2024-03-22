library(tidyverse)
library(furrr)
library(coloc)

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

langefeld_stats <- 
    "./data/langefeld_clean.tsv" |>
    read_tsv() |>
    filter(!is.na(se))

#eqtl_err <- 
#    list.files("./results", pattern = "QTD\\d+\\.err", full.names = TRUE) |>
#    {function(x) setNames(x, str_extract(x, "QTD\\d+"))}() |>
#    map(read_lines)

plan(multisession, workers = availableCores())

eqtl_stats <- 
    list.files("./results", pattern = "QTD\\d+\\.tsv", full.names = TRUE) |>
    {function(x) setNames(x, str_extract(x, "QTD\\d+"))}() |>
    future_map_dfr(~read_tsv(.) |> 
		   select(chr = chromosome, 
			  molecular_trait_id, 
			  gene_id, 
			  pos = position, 
			  rsid, ref, alt, an, beta, maf, se, pvalue),
		   .id = "dataset")

# Match by position and alleles.
# I found that matching by RsID is tricky due to all sorts of problems:
# RsIDs change through time. Same position and alleles are associated with different RsIDs within and across datasets.
# For delins, positions can vary across datasets and dbSNP.
# When there are multiple RsIDs at the same position (different variant species), sometimes RsIDs naming is wrong in some dataset.

# Alleles match in order
min_df <- 
    eqtl_stats |>
    group_by(dataset, gene_id, molecular_trait_id, chr, pos, ref, alt, an, beta, maf, se, pvalue) |>
    summarise(rsid = paste(rsid, collapse = "/")) |>
    ungroup() |>
    inner_join(langefeld_stats,
	       join_by(pos, ref == alleleA, alt == alleleB),
	       suffix = c("_eqtl", "_gwas")) |>
    unite("var_id", c(chr, pos, ref, alt), sep = "-") |>
    select(dataset, gene_id, molecular_trait_id, var_id, 
	   beta_gwas, se_gwas, beta_eqtl, se_eqtl, an_eqtl = an, maf_eqtl = maf)

# Run coloc
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

coloc_res <- min_df |>
    unite("subset_id", c(dataset, gene_id, molecular_trait_id), sep = ",") |> 
    {function(x) split(x, x$subset_id)}() |>
    map(run_coloc)

coloc_results_summary <- map(coloc_res, "summary") |>
    map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
    pivot_wider(names_from = stat, values_from = value) |>
    select(id, nsnps, pp0 = PP.H0.abf, pp1 = PP.H1.abf, pp2 = PP.H2.abf, pp3 = PP.H3.abf, pp4 = PP.H4.abf) |>
    separate(id, c("dataset_id", "gene_id", "molecular_trait_id"), sep = ",")

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

# eQTL catalogue metadata
eqtl_meta <- 
    read_tsv("../data/eqtl_datasets.tsv") |>
    select(-study_id, -tissue_id)

out <- 
    left_join(coloc_results_summary, gene_info, join_by(gene_id)) |>
    left_join(eqtl_meta, join_by(dataset_id)) |>
    select(dataset_id, study = study_label, sample_group, tissue = tissue_label,
	   condition = condition_label, level = quant_method,
	   gene_id, gene_name, molecular_trait_id, everything())

write_tsv(out, "./results/coloc_results.tsv")


## Wang et al
wang <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Wang2021/ASN/harmonized/33536424-GCST90011866-EFO_0002690.h.tsv.gz" |>
    data.table::fread() |>
    filter(!is.na(beta)) |>
    filter(hm_chrom == 7) |> 
    select(chrom = hm_chrom, pos = hm_pos, rsid = hm_rsid,  
	   other_allele = hm_other_allele, effect_allele = hm_effect_allele, 
	   p = p_value, beta, se = standard_error) |>
    arrange(pos) |>
    as_tibble()

min_df_wang <- 
    inner_join(eqtl_stats, wang,
	       join_by(pos, rsid, ref == other_allele, alt == effect_allele),
	       suffix = c("_eqtl", "_gwas")) |>
    select(dataset, gene_id, molecular_trait_id, var_id = rsid, 
	   beta_gwas, se_gwas, beta_eqtl, se_eqtl, an_eqtl = an, maf_eqtl = maf)

coloc_res_wang <- min_df_wang |>
    unite("subset_id", c(dataset, gene_id, molecular_trait_id), sep = ",") |> 
    {function(x) split(x, x$subset_id)}() |>
    map(run_coloc)

coloc_results_summary_wang <- map(coloc_res_wang, "summary") |>
    map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
    pivot_wider(names_from = stat, values_from = value) |>
    select(id, nsnps, pp0 = PP.H0.abf, pp1 = PP.H1.abf, pp2 = PP.H2.abf, pp3 = PP.H3.abf, pp4 = PP.H4.abf) |>
    separate(id, c("dataset_id", "gene_id", "molecular_trait_id"), sep = ",")



