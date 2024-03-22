library(tidyverse)

datasets <- read_tsv("./data/eqtl_datasets.tsv")

temp_dir <- system("echo $TEMP_WORK", intern = TRUE) |>
    file.path("eqtl_catalogue")

query_df <- list.files(temp_dir, full.names = TRUE) |>
    {function(x) setNames(x, sub("^([^_]+).+$", "\\1", basename(x)))}() |>
    map_dfr(read_tsv, .id = "dataset_id")

write_tsv(query_df, "./data/eqtl_catalogue_associations.tsv.gz")


#### test
query_df <- read_tsv("./data/eqtl_catalogue_associations.tsv.gz")

eqtl_vars <- query_df |>
    distinct(chr = chromosome, pos = position, rsid, type, ref, alt)

bentham <- read_tsv("./data/summstats_bentham.tsv") |>
    select(chr, pos, rsid, effect_allele, other_allele, beta, se, p = p_value) |>
    mutate_at(vars(effect_allele, other_allele), toupper)

langefeld <- read_tsv("./data/summstats_langfeld.tsv") |>
    mutate(chr = 1) |>
    select(chr, pos, rsid, effect_allele = alleleA, other_allele = alleleB, beta, se, p = p_lrt) |>
    mutate_at(vars(effect_allele, other_allele), toupper)

wang_eas <- read_tsv("./data/summstats_wang_eas.tsv") |>
    select(chr, pos, rsid, effect_allele, other_allele, beta, se, p = p_value) |>
    mutate_at(vars(effect_allele, other_allele), toupper)
    
wang_eur <- read_tsv("./data/summstats_wang_eur.tsv") |>
    select(chr, pos, rsid, effect_allele = allele1, other_allele = allele2, z = zscore, p = p_value) |>
    mutate_at(vars(effect_allele, other_allele), toupper)

wang_meta <- read_tsv("./data/summstats_wang_meta.tsv") |>
    select(chr, pos, rsid, effect_allele, other_allele, beta, se, p = p_value) |>
    mutate_at(vars(effect_allele, other_allele), toupper)


complement <- function(string) {
    base_compl <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    bases <- unlist(strsplit(string, ""))
    paste(base_compl[bases], collapse = "")
}
    
gwas_stats_list <- 
    list("Bentham" = bentham,
	 "Langefeld" = langefeld,
	 "Wang-EUR" = wang_eur,
	 "Wang-EAS" = wang_eas,
	 "Wang-Meta" = wang_meta)

gwas_stats <- gwas_stats_list |>
    map_dfr(~select(., chr:other_allele), .id = "gwas")


# Filter for alleles that match the eQTL data for both alleles
# Also try to flip alleles to see if they match
# This will not ensure correct effect signs because strand ambiguity and reverse complement alleles across genome assemblies
gwas_stats_match <- gwas_stats |>
    inner_join(eqtl_vars, join_by(chr, rsid), relationship = "many-to-many") |>
    mutate(eff_a_comp = map_chr(effect_allele, complement),
	   oth_a_comp = map_chr(other_allele, complement)) |>
    filter((effect_allele == ref & other_allele == alt) |
	   (other_allele == ref & effect_allele == alt) |
	   (eff_a_comp == ref & oth_a_comp == alt) |
	   (oth_a_comp == ref & eff_a_comp == alt)) |>
    select(gwas, chr, pos = pos.x, rsid, effect_allele, other_allele)

# need to join in a list because different GWAS have different stats
tmp <- gwas_stats_match |>
    {function(x) split(x, x$gwas)}() |>
    {function(x) x[names(gwas_stats_list)]}() |>
    map2(gwas_stats_list, ~inner_join(.x, .y))





