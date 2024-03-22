library(tidyverse)
library(susieR)

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

vcf <- 
    "./data/chr7:49766324-50738073.vcf.gz" |>
    read_tsv(comment = "##") |>
    rename("chr" = `#CHROM`) |>
    mutate(chr = str_remove(chr, "chr")) |>
    unite("ID", c(chr, POS, REF, ALT), sep = "_", remove = FALSE) |>
    select(chr, POS, ID, everything())

plink <- 
    "./data/chr7:49766324-50738073.ld" |>
    data.table::fread(header = FALSE) |>
    as_tibble()

colnames(plink) <- vcf$ID

valid_snps <- sapply(plink, function(x) !any(is.na(x)))

plink_filt <- plink[valid_snps, valid_snps]
ldmat <- data.matrix(plink_filt)
rownames(ldmat) <- colnames(ldmat)
vcf <- filter(vcf, valid_snps)

langefeld_sample_size <- 6748 + 11516

langefeld_stats <- 
    "./data/langefeld_clean.tsv" |>
    read_tsv() |>
    mutate(chr = "7") |>
    select(chr, rsid, pos, ref, alt, beta, se, p = p_lrt) |>
    unite("ID", c(chr, pos, ref, alt), sep = "_", remove = FALSE) |>
    filter(ID %in% vcf$ID) |>
    drop_na(se) |>
    mutate(z = beta/se,
	   ID = factor(ID, levels = vcf$ID)) |>
    arrange(ID)

ldmat <- ldmat[langefeld_stats$ID, langefeld_stats$ID]

condz <- kriging_rss(langefeld_stats$z, ldmat, n = langefeld_sample_size, r_tol = 1e-04)
png("./plots/diagnostic_susie.png")
condz$plot
dev.off()

snp_rm <- 
    as_tibble(condz$conditional_dist) |>
    rowid_to_column() |>
    arrange(desc(abs(z_std_diff))) |>
    filter(z_std_diff >= 6)

colnames(ldmat)[snp_rm$rowid]


# if it does not converge, remove problematic SNPs and repeat until convergence
iter <- 0L
fit <- list()
fit$converged <- FALSE
while ( !fit$converged & iter <= 20 ) {

    if ( iter > 0 ) { 
    
	condz <- kriging_rss(summ_stats$z, ldmat, n = sample_size, r_tol = 1e-04)

	snp_rm <- as_tibble(condz$conditional_dist) |>
	    rowid_to_column() |>
	    arrange(desc(abs(z_std_diff))) |>
	    slice(1) |>
	    pull(rowid)

	ldmat <- ldmat[-snp_rm, -snp_rm]
	summ_stats <- slice(summ_stats, -snp_rm)
    }

    fit <- susie_rss(summ_stats$z, 
		     ldmat, 
		     n = sample_size, 
		     L = 10, 
		     coverage = 0.9, 
		     r_tol = 1e-05)

    iter <- iter + 1L
}


make_pip_df <- function(fit_obj) {

    if (! is.null(fit_obj$sets$coverage) ) {

	coverage_df <- 
	    tibble(coverage = scales::percent(round(fit_obj$sets$coverage, 3)),
		   cs = paste0("L", fit_obj$sets$cs_index))

	cs_df <- enframe(fit_obj$sets$cs, "cs", "rowid") |>
	    unnest(cols = rowid) |>
	    left_join(coverage_df, by = "cs")

	pip_df <- enframe(fit_obj$pip, "ID", "pip") |>
	    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
	    rowid_to_column() |>
	    left_join(cs_df) |>
	    mutate(cs = ifelse(!is.na(cs), paste0(cs, " (", coverage, ")"), NA))

    } else {
	
	pip_df <- 
	    enframe(fit_obj$pip, "ID", "pip") |>
	    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
	    mutate(cs = NA)
    }

    return(pip_df)
}



# Analysis
plan(multisession, workers = availableCores())

susie_results_langefeld <- langefeld_regions$locus |>
    setNames(langefeld_regions$locus) |>
    future_map(run_susie, 
	       region_df = langefeld_regions, 
	       summ_stats_df = langefeld_summ_stats_all, 
	       sample_size = langefeld_sample_size)

write_rds(susie_results_langefeld, "./susie_results_langefeld.rds")

langefeld_pip_df <- map_dfr(susie_results_langefeld, make_pip_df, .id = "locus")

write_tsv(langefeld_pip_df, "./susie_results_langefeld.tsv")

