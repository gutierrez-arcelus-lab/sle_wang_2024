library(tidyverse)
library(susieR)

# Process VCF and LD matrix
Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)
vcf <- "./data/chr1_eur.vcf.gz" |>
    read_tsv(comment = "##") |>
    unite("ID", c("ID", "REF", "ALT"), sep = "_")

plink <- "./data/chr1_eur.ld" |>
    read_delim(delim = " ", col_names = FALSE)

if ( (ncol(plink) %% nrow(plink) == 1) && all(is.na(plink[, ncol(plink)])) ) {

    plink <- plink[, -ncol(plink)]
}

colnames(plink) <- vcf$ID
plink <- plink |> 
    add_column("rsid" = vcf$ID, .before = 1)

valid_snps <- which(sapply(plink, function(x) !all(is.na(x))))

plink_filt <- plink[valid_snps, valid_snps]
####

# GWAS data
bentham_n <- 4036 + 6959

bentham_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/opengwas/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##")

snp_window <- bentham_stats |> 
    filter(ID == "rs2297550") |>
    mutate(start = POS - 5e5L, stop = POS + 5e5L) |>
    select(`#CHROM`, start, stop)

bentham_stats_window <- bentham_stats |>
    inner_join(snp_window, join_by(`#CHROM`, between(POS, start, stop))) |>
    separate_rows(FORMAT:`EBI-a-GCST003156`, sep = ":") |>
    select(chr = `#CHROM`, pos = POS, rsid = ID, REF, ALT, FORMAT, stats = `EBI-a-GCST003156`) |>
    pivot_wider(names_from = FORMAT, values_from = stats) |>
    mutate(ID = paste(ID, REF, ALT, sep = "_")) |>
    select(chr, pos, ID, beta = ES, se = SE) |>
    mutate_at(vars(beta, se), as.numeric)

bentham_final_stats <- bentham_stats_window |>
    filter(ID %in% plink$rsid) |>
    mutate(z = beta/se)

plink_final <- plink_filt |>
    filter(rsid %in% bentham_final_stats$ID) |>
    select(rsid, bentham_final_stats$ID) |>
    select(-rsid)

ldmat <- data.matrix(plink_final)
rownames(ldmat) <- colnames(ldmat)

condz <- kriging_rss(bentham_final_stats$z, ldmat, n = bentham_n, r_tol = 1e-04)

dev_df <- 
    condz$conditional_dist |>
    as_tibble() |>
    add_column(rsid = colnames(ldmat), .before = 1) |>
    filter(logLR > 2 & abs(z) > 2) 

p <- condz$conditional_dist |>
    as_tibble() |>
    add_column(rsid = colnames(ldmat), .before = 1) |>
    separate(rsid, c("ID", "REF", "ALT"), sep = "_") |>
    left_join(select(bentham_stats, ID, INFO), join_by(ID)) |>
    ggplot(aes(x = condmean, y = z, color = INFO)) +
	geom_abline() +
	geom_point() +
	scale_color_manual(values = c("." = "black", "ReverseComplementedAlleles" = "red")) +
	ggrepel::geom_text_repel(data = dev_df, 
				 aes(x = condmean, y = z, label = rsid), 
				 size = 2.5, inherit.aes = FALSE) +
	theme(legend.position = "top")


ggsave("./plots/susie_diagnostic.png", p)

fit <- susie_rss(bentham_final_stats$z, 
		 ldmat, 
		 n = bentham_n, 
		 L = 3, 
		 coverage = 0.5, 
		 r_tol = 1e-05)

fit$converged
# if it does not converge, remove problematic SNPs and repeat until convergence

coverage_df <- 
    tibble(coverage = scales::percent(round(fit$sets$coverage, 3)),
	   cs = paste0("L", fit$sets$cs_index))

#cs_df <- enframe(fit_obj$sets$cs, "cs", "rowid") |>
#    unnest(cols = rowid) |>
#    left_join(coverage_df, by = "cs")
#
#pip_df <- enframe(fit_obj$pip, "ID", "pip") |>
#    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
#    rowid_to_column() |>
#    left_join(cs_df) |>
#    mutate(cs = ifelse(!is.na(cs), paste0(cs, " (", coverage, ")"), NA))


enframe(fit$pip, "ID", "pip") |>
    separate(ID, c("rsid", "ref", "alt"), sep = "_") |>
    mutate(cs = NA) |>
    arrange(desc(pip))
