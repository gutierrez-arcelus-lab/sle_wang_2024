library(tidyverse)
library(rvest)
library(patchwork)
library(httr)
library(glue)
library(jsonlite)



# Summary stats
bentham <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690-build37.f.tsv.gz" |>
    read_tsv(col_types = c(chromosome = "c")) |>
    janitor::clean_names() |>
    select(chr = chromosome, pos = base_pair_location, rsid = variant_id,
	   effect_allele, other_allele, beta, se = standard_error, p_value) |>
    filter(chr == 1, between(pos, 200e6, 210e6))

wang_eur <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Wang2021/EUR/Meta_Results.txt" |>
    read_tsv(col_types = c(CHR = "c")) |>
    janitor::clean_names() |>
    select(chr, pos = bp, rsid = snp, allele1 = a1lele1, allele2, zscore, p_value) |>
    filter(chr == 1, between(pos, 200e6, 210e6))

wang_asn <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Wang2021/ASN/GCST90011866_buildGRCh37.tsv" |>
    read_tsv(col_types = c(chromosome = "c")) |>
    select(chr = chromosome, pos = base_pair_location, rsid = variant_id,
	   effect_allele, other_allele, beta, se = standard_error, p_value) |>
    filter(chr == 1, between(pos, 200e6, 210e6))

wang_meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Wang2021/meta/GCST011096_buildGRCh37.tsv" |>
    read_tsv(col_types = c(chromosome = "c")) |>
    select(chr = chromosome, pos = base_pair_location, rsid = variant_id,
	   effect_allele, other_allele, beta, se = standard_error, p_value) |>
    filter(chr == 1, between(pos, 200e6, 210e6))

langefeld_chr1 <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr1.out" |>
    read_delim(comment = "#", delim = " ") |>
    select(rsid, 
	   pos = position, 
	   alleleA, 
	   alleleB,
	   beta = `frequentist_add_beta_1:add/sle=1`, 
	   se = frequentist_add_se_1, 
	   p_lrt = frequentist_add_lrt_pvalue,
	   p_wald = frequentist_add_wald_pvalue_1) |>
    drop_na(beta) |> 
    filter(grepl("^rs", rsid)) |>
    extract(rsid, "rsid", "(rs\\d+)") |>
    filter(between(pos, 200e6, 210e6))

japanese <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Japanese_SLE_20210508/METAL_Auto_X_with_GC.meta" |>
    read_delim(delim = " ") |>
    filter(grepl("^1:", MARKER)) |>
    extract(MARKER, c("chr", "pos"), "([^:]+):(\\d+)") |>
    mutate(pos = as.numeric(pos),
	   rsid = NA) |>
    select(chr, pos, rsid, effect_allele = ALT, other_allele = REF, beta = ALT_BETA, se = SEBETA, p_value = PVALUE)

# Save summ stats
write_tsv(bentham, "./data/summstats_bentham.tsv")
write_tsv(wang_eur, "./data/summstats_wang_eur.tsv")
write_tsv(wang_asn, "./data/summstats_wang_eas.tsv")
write_tsv(wang_meta, "./data/summstats_wang_meta.tsv")
write_tsv(langefeld_chr1, "./data/summstats_langfeld.tsv")
write_tsv(japanese, "./data/summstats_japanese.tsv")

# selected window to plot
snp <- 
    filter(bentham, rsid == "rs2297550") |>
    mutate(start = pos - 5e5, stop = pos + 5e5) |>
    select(chr, start, stop)

gwas_stats <- 
    bind_rows(
	      "Bentham et al. (EUR)" = 
		  inner_join(bentham, snp, join_by(chr, between(pos, start, stop))) |>
		    select(rsid, pos, effect_allele, other_allele, p = p_value),
	      "Langefeld et al. (EUR)" = 
		  inner_join(langefeld_chr1, snp, join_by(between(pos, start, stop))) |>
		    select(rsid, pos, effect_allele = alleleA, other_allele = alleleB, p = p_wald),
	      "Wang et al. (EUR)" = 
		  inner_join(wang_eur, snp, join_by(chr, between(pos, start, stop))) |>
		    select(rsid, pos, effect_allele = allele1, other_allele = allele2, p = p_value),
	      "Wang et al. (EAS)" = 
		  inner_join(wang_asn, snp, join_by(chr, between(pos, start, stop))) |>
		    select(rsid, pos, effect_allele, other_allele, p = p_value),
	      "Wang et al. (EUR+EAS meta-analysis)" = 
		  inner_join(wang_meta, snp, join_by(chr, between(pos, start, stop))) |>
		    select(rsid, pos, effect_allele, other_allele, p = p_value),
	      "Japanese" = 
		  inner_join(japanese, snp, join_by(chr, between(pos, start, stop))) |>
		    select(rsid, pos, effect_allele, other_allele, p = p_value),
	      .id = "GWAS") |>
    mutate(GWAS = fct_inorder(GWAS)) |>
    drop_na(p) |>
    mutate_at(vars(effect_allele, other_allele), toupper)

# Save variants in GWAS
gwas_stats |>
    arrange(pos) |>
    distinct(rsid) |>
    filter(grepl("^rs", rsid)) |>
    pull(rsid) |>
    write_lines("./data/gwas_rsid.txt")

# Plot regions
gwas_plot <-     
    ggplot(data = gwas_stats, 
	   aes(x = pos, y = -log10(p))) +
    geom_vline(xintercept = filter(bentham, rsid == "rs2297550") |> pull(pos),
	       linewidth = .2, linetype = 2) +
    geom_vline(xintercept = 
	       filter(bentham, rsid == "rs2297550") |> 
	       pull(pos) |> 
	       {function(x) x + c(-2.5e5, 2.5e5)}(),
	       linewidth = .2, linetype = 2, color = "blue") +
    geom_point(size = .75, alpha = .5) +
    geom_point(data = filter(gwas_stats, pos == 206643772),
	       aes(x = pos, y = -log10(p)),
	       size = 1, shape = 1, color = "tomato3") +
    geom_label(data = filter(gwas_stats, pos == 206643772) |> fill(rsid),
	       aes(x = pos, y = -log10(p), label = rsid),
	       alpha = .75, label.padding = unit(0.05, "lines"),
	       size = 2.5, color = "tomato3", hjust = 0, vjust = 0) +
    facet_wrap(~GWAS, ncol = 1) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 8),
	  axis.title.y = element_text(size = 9),
	  axis.ticks = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.y = element_line(linewidth = .25),
	  panel.grid.minor.y = element_blank(),
	  strip.text = element_text(face = "bold"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL,
	 y = "log10 P-value") +
    coord_cartesian(clip = "off")


# Gene tracks
tracks <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gene_tracks_gencodev19.tsv" |>
    read_tsv()
  
# fix repeated gene names
track_names <- tracks |>
    distinct(gene_id, gene_name) |>
    add_count(gene_name) |>
    group_by(gene_name) |>
    mutate(i = seq_len(n()),
	   gene_name = ifelse(n > 1, paste(gene_name, i, sep = "-"), gene_name)) |>
    ungroup() |>
    select(gene_id, gene_name) |>
    arrange(gene_name)

track_df <- 
    tracks |>
    filter(chr == paste0("chr", snp$chr),
	   between(start, snp$start, snp$stop) | between(end, snp$start, snp$stop)) |> 
    select(chr, gene_id, start, end, strand, model = feature) |>
    left_join(track_names) |>
    mutate(start = ifelse(start < min(gwas_stats$pos), min(gwas_stats$pos), start),
	   end = ifelse(end > max(gwas_stats$pos), max(gwas_stats$pos), end))

n_gene_groups <- 4

gene_ys <- track_df |>
    group_by(gene_name) |>
    summarise(start = min(start)) |>
    ungroup() |>
    arrange(start) |>
    mutate(g = ntile(start, n_gene_groups)) |>
    group_by(g) |>
    mutate(y = seq_len(n()) + 1L) |>
    ungroup() |>
    select(gene_name, y)

track_df <- left_join(track_df, gene_ys)

gene_labels <- track_df |>
    group_by(gene_name, y) |>
    summarise(s = min(start)) |>
    ungroup()

arrows_df <- track_df |>
    group_by(gene_name) |>
    mutate(i = ifelse(strand == "+", max(end), min(start))) |>
    ungroup() |>
    distinct(gene_name, start = i, y, strand) |>
    mutate(end = ifelse(strand == "+", start + 1e2, start - 1e2))

gene_plot <- 
    ggplot(track_df) +
    geom_segment(aes(x = start, xend = end, y = y, yend = y, linewidth = model),
		 color = "midnightblue") +
    geom_segment(data = arrows_df,
		 aes(x = start, xend = end, y = y, yend = y),
		 arrow = arrow(length = unit(.15, "cm")),
		 color = "midnightblue") +
    geom_text(data = gene_labels, 
	      aes(x = s, y = y, label = gene_name),
	      hjust = 1.15, size = 2, color = "grey40", fontface = "italic") +
    scale_x_continuous(limits = range(gwas_stats$pos),
				labels = function(x) x/1e6L) +
    scale_y_continuous(breaks = 0:max(gene_ys$y)) +
    scale_linewidth_manual(values = c("exon" = 4, "intron" = 1)) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
	  axis.ticks.y = element_blank(),
	  axis.line.x = element_line(color = "black", linewidth = .25),
	  axis.text.x = element_text(size = 8),
	  axis.title.x = element_text(size = 9),
	  legend.position = "none",
	  panel.border = element_rect(color = NA),
	  plot.background = element_rect(fill = "white", color = "white"),
	  plot.margin = margin(r = 1, l = 1, unit = "cm"),
	  panel.grid = element_blank()) +
    labs(x = "Position on chr1 hg19 (Mb)", y = NULL) +
    coord_cartesian(clip = "off")

ggsave("./plots/gwas.png", 
       gwas_plot + gene_plot + plot_layout(heights = c(1, .2)),
       height = 7, width = 6, dpi = 600)


# 1000 Genomes
kgp_pops <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv" |>
    read_tsv() |>
    janitor::clean_names() |>
    filter(population_description != "Total")

kgp_samples <- 
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##") |>
    select(sample_name = SAMPLE_NAME, population_code = POPULATION) |>
    distinct() |>
    left_join(kgp_pops, by = "population_code") |>
    select(sample_name, population_code, super_population)

kgp_eur_eas <- kgp_samples |>
    filter(super_population %in% c("EUR", "EAS")) |>
    arrange(sample_name)

kgp_eur_eas |>
    pull(sample_name) |>
    write_lines("./data/kgp_eur_eas.txt")

kgp_eur_eas |>
    filter(super_population == "EUR") |>
    pull(sample_name) |>
    write_lines("./data/kgp_eur.txt")

kgp_eur_eas |>
    filter(super_population == "EAS") |>
    pull(sample_name) |>
    write_lines("./data/kgp_eas.txt")



# eQTL Catalogue
#All datasets will be pulled if this parameter is bigger than the actual number of datasets
max_pulled_rows <- 1000 

URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")

# Make a request
r <- GET(URL, accept_json())

# Check status
status_code(r)

# Extract content and convert content to dataframe
datasets <- 
    content(r, "text", encoding = "UTF-8") |>
    fromJSON() |>
    as_tibble()

#dir.create("data")

write_tsv(datasets, "./data/eqtl_datasets.tsv")

