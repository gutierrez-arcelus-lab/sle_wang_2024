library(tidyverse)

coloc_res <- 
    read_tsv("./results/coloc_results.tsv") |>
    mutate(study = ifelse(grepl("_\\d+$", study), paste0(sub("_", " (", study), ")"), study),
	   condition = sub("_", " ", condition),
	   study_label = ifelse(condition == "naive", 
				paste(study, tissue), 
				paste(study, tissue, condition)),
	   level = recode(level, "ge" = "gene"),
	   level = factor(level, levels = c("gene", "exon", "microarray", "aptamer", "tx", "txrev", "leafcutter"))) |>
    select(study_label, level, gene_id, gene_name, molecular_trait_id, pp4)

coloc_filt <- 
    coloc_res |> 
    group_by(study_label, level, gene_id, gene_name) |>
    slice_max(pp4) |>
    group_by(study_label) |>
    filter(any(pp4 > .5)) |>
    group_by(gene_id) |>
    filter(any(pp4 > .1)) |>
    ungroup()



npg_colors <- c(ggsci::pal_npg()(10)[c(2, 6, 5, 8)], "black")

pl <- 
    ggplot(coloc_filt, aes(x = pp4, y = study_label)) +
    geom_point(aes(color = gene_name)) +
    scale_color_manual(values = npg_colors)  +
    scale_x_continuous(breaks = c(0, 1),
		       labels = c("0", "1")) +
    facet_wrap(~level, nrow = 1) +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "PP4", y = NULL, color = "Gene: ")

ggsave("./plots/coloc.png", pl, width = 8, height = 10)


#
library(AnnotationHub)
library(locuszoomr)

#data(SLE_gwas_sub)

risk_pos <- 50266267 

ah <- AnnotationHub()
#query(ah, "EnsDb.Hsapiens.v105")

ens_data <- ah[["AH98047"]]

ld_token <- read_lines("../ldlink_api_token.txt")

langefeld <- 
    "./data/langefeld_clean.tsv" |>
    data.table::fread() |>
    dplyr::mutate(chrom = 7) |>
    dplyr::select(chrom, pos, rsid, other_allele = alleleA, effect_allele = alleleB,
		  p = p_wald, beta, se)

bentham <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    data.table::fread() |>
    dplyr::filter(!is.na(beta)) |>
    dplyr::filter(hm_chrom == 7) |> 
    dplyr::select(chrom = hm_chrom, pos = hm_pos, rsid = hm_rsid,  
	   other_allele = hm_other_allele, effect_allele = hm_effect_allele, 
	   p = p_value, beta, se = standard_error)

loc <- 
    locus(data = bentham, 
	  gene = 'IKZF1',
	  flank = 1e5,
	  ens_db = ens_data)

loc <- link_LD(loc, token = ld_token)

loc_ggplot <- locus_ggplot(loc)

