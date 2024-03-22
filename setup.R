library(tidyverse)
library(glue)

gene_annots <-
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "##", col_names = FALSE) |>
    filter(X3 == "gene")

# Extract IKZF1 summary stats from eQTL Catalogue
cat_dir <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue"

ge_files <-
    file.path(cat_dir, "ge_file_list.txt") |>
    read_lines() |>
    basename()

ikzf1_tss <- 
    gene_annots |>
    filter(X1 == "chr7") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    filter(gene_name == "IKZF1") |>
    pull(X4)

ikzf1_window <- 
    paste(ikzf1_tss + c(-1e6, 1e6), collapse = "-") |>
    {function(x) paste0("7:", x)}()

out_dir <- "./data/eQTLcatalogue/sumstats/ge" 
if (!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ge_files_full <- file.path(cat_dir, "ge_all", ge_files)
dataset_ids <- sub("^(QTD\\d+)\\..+$", "\\1", ge_files)

cmd <- 
    glue("tabix {ge_files_full} {ikzf1_window} | bgzip > {out_dir}/{dataset_ids}.ikzf1.tsv.gz")

walk(cmd, system)
