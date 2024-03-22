library(tidyverse)
library(ggsci)

coloc_df <- read_tsv("./data/coloc_results.tsv")

plot_df <- coloc_df |>
    mutate(tissue_label = case_when(condition == "naive" ~ tissue,
				    TRUE ~ paste(tissue, condition))) |>
    group_by(study, sample_group, tissue_label, gene_id, level) |>
    slice_max(h4) |>
    ungroup() |>
    mutate(study_id = sprintf("%s (%s)", study, tissue_label),
	   study_id = gsub("_", " ", study_id)) |>
    select(study_id, level, gene_id, gene_name, pp4 = h4) |>
    group_by(gene_id) |>
    filter(any(pp4 > .2)) |>
    group_by(study_id) |>
    filter(any(pp4 > .2)) |>
    ungroup() |>
    mutate(level = factor(level, levels = c("ge", "exon", "tx", "txrev", "leafcutter"))) |>
    arrange(desc(pp4)) |>
    mutate(study_id = fct_inorder(study_id),
	   gene_name = fct_relevel(gene_name, "IKBKE", after = 0))


p <- ggplot(plot_df, aes(x = pp4, y = study_id)) +
    geom_point(aes(color = gene_name), size = 3) +
    facet_wrap(~level, nrow = 1) +
    scale_x_continuous(limits = c(0, 1),
		       breaks = c(0, 1),
		       labels = c("0", "1")) +
    scale_color_npg() +
    theme(panel.background = element_rect(fill = "grey96")) +
    labs(x = "Posterior probability of colocalization", 
	 y = NULL,
	 color = "Gene")

ggsave("./plots/colocs.png", p, width = 10, height = 7, dpi = 600)
