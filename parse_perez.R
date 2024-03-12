library(tidyverse)

eqtls <- 
    "./data/perez_et_al/science.abf1970_table_s6.txt" |>
    read_tsv(skip = 1) |>
    select(-STD_FE) |>
    extract(RSID, c("chr", "pos", "gene"), "([^:]+):([^_]+)_(.+)", convert = TRUE)


eqtls |> filter(chr == 1, pos > 206e6)

eqtls |>
    filter(chr == 1, pos == 206470429)

eqtls |>
    filter(chr == 1, pos == 206643772)

"./data/perez_et_al/science.abf1970_table_s10.txt" |>
    read_tsv(skip = 1) |>
    filter(gene == "IKZF1") |>
    arrange(desc(H4))
