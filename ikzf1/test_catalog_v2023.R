library(tidyverse)
library(httr)
library(glue)
library(jsonlite)

max_pulled_rows <- 1000
URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")

r <- GET(URL, accept_json())
status_code(r)

cont <- content(r, "text", encoding = "UTF-8")
datasets <- 
    fromJSON(cont) |>
    as_tibble()

write_tsv(datasets, "catalogue_datasets_Mar2024.tsv")

datasets |>
    filter(grepl("macrophage|monocyte|neutrophil|blood|T-cell|T cell|Treg|Th\\d+|Tfh|B cell|LCL|NK cell", tissue_label)) |>
    write_tsv("catalogue_datasets_Mar2024_immune.tsv")

datasets |> 
    filter(grepl("Schmiedel", study_label)) |> 
    filter(tissue_label == "CD8+ T cell", condition_label == "naive") |>
    print(width = Inf)

read_tsv("results/QTD000492.err") |> 
    distinct(gene_id)

dat <- read_tsv("./QTD000492.cc.tsv.gz")

dat |> 
    filter(chromosome == 7, position > 46e6) |> arrange(position)
    distinct(molecular_trait_id) |>
    extract(molecular_trait_id, "gene_id", "(ENSG\\d+)\\..+") |>
    distinct() |> filter(gene_id == "ENSG00000185811")
