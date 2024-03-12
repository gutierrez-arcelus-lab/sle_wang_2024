library(tidyverse)
library(coloc)

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

langefeld_stats <- 
    "./data/langefeld_clean.tsv" |>
    read_tsv() |>
    filter(!is.na(se)) |>
    mutate(chr = "7") |>
    select(chr, pos, rsid, allele_a = alleleA, allele_b = alleleB, beta, se)

bentham_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    read_tsv() |>
    drop_na(beta) |>
    filter(hm_chrom == 7, 
	   between(hm_pos, min(langefeld_stats$pos), max(langefeld_stats$pos))) |>
    select(chr = hm_chrom, rsid = hm_rsid, pos = hm_pos,
	   allele_a = hm_other_allele, allele_b = hm_effect_allele, 
	   beta, se = standard_error)

min_df <- 
    inner_join(langefeld_stats, bentham_stats,
    join_by(chr, pos, rsid, allele_a, allele_b),
    suffix = c("_langef", "_bentham")) |>
    filter(beta_langef != 0, beta_bentham != 0)

# Run coloc
gwas_1 <- 
    list(beta = min_df$beta_langef,
	 varbeta = min_df$se_langef^2,
	 type = "cc",
	 snp = min_df$rsid)

gwas_2 <- 
    list(beta = min_df$beta_bentham,
	 varbeta = min_df$se_bentham^2,
	 type = "cc",
	 snp = min_df$rsid)

coloc_res <- coloc.abf(dataset1 = gwas_1, 
		       dataset2 = gwas_2)

coloc_res
