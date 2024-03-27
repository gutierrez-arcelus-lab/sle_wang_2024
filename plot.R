library(tidyverse)
library(extrafont)
library(glue)

"./data/eQTLcatalogue/metadata.tsv" |>
    read_tsv() |>
    distinct(tissue_label) |>
    pull(tissue_label)

immune <- 
    c("macrophage", "monocyte", "neutrophil", "CD4+ T cell", "Treg memory", "LCL", 
      "CD8+ T cell", "platelet", "B cell", "T cell", "NK cell", "blood", "Tfh cell",
      "Th17 cell", "Th1 cell", "Th2 cell", "Treg naive", "CD16+ monocyte", "synovium", "microglia")

coloc_res <- 
    read_tsv("./results/coloc_results_ikzf1.tsv") |>
    filter(tissue_label %in% immune) |>
    mutate(study_label = sub("(_)(\\d+)$", " (\\2)", study_label),
	   condition_label = sub("_", " ", condition_label),
	   lab = ifelse(condition_label == "naive", 
			paste(study_label, tissue_label), 
			paste(study_label, tissue_label, condition_label)),
	   quant_method = recode(quant_method, "ge" = "gene"),
	   quant_method = factor(quant_method, levels = c("gene", "exon"))) |>
    select(gwas_id, lab, quant_method, gene_id, gene_name, molecular_trait_id, pp4)

coloc_filt <- 
    coloc_res |> 
    group_by(gwas_id, lab, quant_method, gene_id, gene_name) |>
    slice_max(pp4) |>
    group_by(gwas_id, lab, quant_method) |>
    filter(any(pp4 > .5)) |>
    ungroup() |>
    mutate(gene_name = ifelse(pp4 > 0.1, gene_name, "Other")) |>
    select(-gene_id, -molecular_trait_id) |>
    arrange(desc(pp4))

my_colors <- c("firebrick", "Deep Sky Blue", "black", "goldenrod", "grey")
names(my_colors) <- unique(coloc_filt$gene_name)

pl <- 
    ggplot(coloc_filt, aes(x = pp4, y = lab)) +
    geom_vline(xintercept = .8, linetype = 2, linewidth = .25) +
    geom_point(aes(color = gene_name), size = 2) +
    scale_color_manual(values = my_colors)  +
    scale_x_continuous(limits = c(0, 1),
		       breaks = c(0, 1),
		       labels = c("0", "1")) +
    facet_grid(quant_method~gwas_id, scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(axis.text = element_text(family = "Arial"),
	  axis.title = element_text(family = "Arial"),
	  legend.text = element_text(family = "Arial"),
	  legend.title = element_text(family = "Arial"),
	  strip.text.x = element_text(family = "Arial"),
	  strip.text.y = element_text(angle = 0, family = "Arial")) +
    labs(x = "Posterior probability of colocalization", 
	 y = NULL, 
	 color = "Gene: ")

ggsave("./plots/coloc.png", pl, width = 6.6, height = 3)


#
library(AnnotationHub)
library(locuszoomr)
library(patchwork)

ah <- AnnotationHub()
#query(ah, "EnsDb.Hsapiens.v105")

ens_data <- ah[["AH98047"]]

risk_var <- "rs4917014"

ld_vcf <- 
    "./data/chr7:49305863-51305863.vcf.gz" |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    unite("ID", c(ID, REF, ALT), sep = "-")

ld_plink <- 
    "./data/chr7:49305863-51305863.ld" |>
    data.table::fread() |>
    as_tibble() |>
    setNames(ld_vcf$ID) |>
    add_column(snp_id = ld_vcf$ID, .before = 1)

ld_risk_var <-
    ld_plink |>
    dplyr::filter(grepl(risk_var, snp_id)) |>
    pivot_longer(-snp_id, names_to = "var_id", values_to = "r2") |>
    dplyr::select(var_id, r2) |>
    separate(var_id, c("rsid", "ref", "alt"), sep = "-", convert = TRUE)

langefeld <- 
    "./data/gwas/langefeld_ikzf1.tsv" |>
    data.table::fread() |>
    dplyr::mutate(chrom = "7") |>
    dplyr::select(chrom, pos, rsid, other_allele = alleleA, effect_allele = alleleB, p)

langefeld_grch38 <- 
    "./data/gwas/langefeld_hg38.bed" |>
    read_tsv(col_names = FALSE) |>
    dplyr::select(chrom = X1, pos = X3, rsid = X4) |>
    dplyr::mutate(chrom = str_remove(chrom, "chr")) |>
    distinct()

langefeld_ld <-
    langefeld |>
    left_join(ld_risk_var, join_by(rsid, other_allele == ref, effect_allele == alt)) |>
    group_by(chrom, pos, rsid) |>
    nest() |>
    ungroup() |>
    inner_join(langefeld_grch38, join_by(chrom, rsid)) |>
    dplyr::select(chrom, pos = pos.y, rsid, data) |>
    unnest(cols = data) |>
    data.table::as.data.table()

loc_langefeld <- 
    locus(data = langefeld_ld, 
	  gene = 'IKZF1',
	  flank = 5e5,
	  LD = "r2",
	  ens_db = ens_data)

loc_langefeld_ggplot <- 
    gg_scatter(loc_langefeld, 
	       labels = "index", 
	       nudge_y = .5, 
	       legend_pos = "topright") +
    labs(title = "Langefeld et al.")

# eQTL Catalogue
dataset_id <- 
    read_tsv("./data/eQTLcatalogue/metadata.tsv") |>
    dplyr::filter(study_label == "GENCORD", tissue_label == "LCL", quant_method == "ge") |>
    pull(dataset_id)

qtl_names <-
    glue("zcat /lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/ge_all/{dataset_id}.all.tsv.gz | head -n1") |>
    system(intern = TRUE) |>
    str_split("\t") |>
    unlist() |>
    str_remove("\\r$")

qtl_data <- 
    glue("./data/eQTLcatalogue/sumstats/{dataset_id}.ikzf1.tsv.gz") |>
    read_tsv(col_names = qtl_names) |>
    dplyr::filter(gene_id == "ENSG00000185811") |>
    dplyr::select(chrom = chromosome, pos = position, rsid, ref, alt, p = pvalue) |>
    dplyr::mutate(chrom = as.character(chrom))

qtl_ld_vcf <- 
    "./data/chr7:49304068-51304068_grch38.vcf.gz" |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    unite("ID", c(`#CHROM`, POS, REF, ALT), sep = "-")

qtl_ld_plink <- 
    "./data/chr7:49304068-51304068_grch38.ld" |>
    data.table::fread() |>
    as_tibble() |>
    setNames(qtl_ld_vcf$ID) |>
    add_column(snp_id = qtl_ld_vcf$ID, .before = 1)

qtl_ld_top_var <-
    qtl_ld_plink |>
    dplyr::filter(grepl("chr7-50266267", snp_id)) |>
    pivot_longer(-snp_id, names_to = "var_id", values_to = "r2") |>
    dplyr::select(var_id, r2) |>
    separate(var_id, c("chrom", "pos", "ref", "alt"), 
             sep = "-", convert = TRUE, remove = FALSE) |>
    dplyr::mutate(chrom = str_remove(chrom, "chr"))

qtl_data_ld <- 
    qtl_data |>
    left_join(qtl_ld_top_var, join_by(chrom, pos, ref, alt)) |>
    dplyr::select(chrom, pos, rsid, ref, alt, p, r2) |>
    data.table::as.data.table()

loc_qtl <- 
    locus(data = qtl_data_ld, 
	  gene = 'IKZF1',
	  flank = 5e5,
	  LD = "r2",
	  ens_db = ens_data)

loc_qtl <- link_recomb(loc_qtl, genome = "hg38")

loc_qtl_ggplot <- 
    gg_scatter(loc_qtl, 
	       labels = "index", 
	       nudge_y = .5, 
	       legend_pos = "topright") +
    labs(title = "GENCORD LCL IKZF1 eQTLs") +
    theme(legend.position = "none")

# Gene tracks
g <- gg_genetracks(loc_langefeld)

# plot
ggsave("./plots/locuszoom.png", 
       loc_langefeld_ggplot / loc_qtl_ggplot / g,
       width = 10, height = 8)


