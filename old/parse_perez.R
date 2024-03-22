library(tidyverse)
library(AnnotationHub)
library(locuszoomr)
library(patchwork)
library("EnsDb.Hsapiens.v75")
library(extrafont)

eqtls <- 
    "./data/perez_et_al/science.abf1970_table_s6.txt" |>
    data.table::fread(skip = "cell") |>
    dplyr::select(-STD_FE) |>
    extract(RSID, c("chrom", "pos", "gene"), "([^:]+):([^_]+)_(.+)", convert = TRUE) |>
    dplyr::filter(gene == "IKBKE") |>
    mutate(chrom = as.character(chrom)) |>
    dplyr::select(cell, chrom, pos, p = PVALUE_FE) |> 
    {function (x) split(x, x$cell)}() |>
    map(~dplyr::select(., -cell))

risk_pos <- 206643772 

ld_vcf <- 
    "./data/chr1:205643772-207643772.vcf.gz" |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    unite("id", c(`#CHROM`, POS, REF, ALT), sep = ":", remove = FALSE) |>
    dplyr::select(chrom = `#CHROM`, pos = POS, rsid = ID, id, ref = REF, alt = ALT) |>
    mutate(chrom = as.character(chrom))

ld_plink <- 
    "./data/chr1:205643772-207643772.ld" |>
    data.table::fread() |>
    as_tibble() |>
    setNames(ld_vcf$id) |>
    add_column(snp_id = ld_vcf$id, .before = 1)

ld_risk_var <-
    ld_plink |>
    dplyr::filter(grepl(paste0("^1:", risk_pos, ":"), snp_id)) |>
    pivot_longer(-snp_id, names_to = "var_id", values_to = "r2") |>
    dplyr::select(var_id, r2) |>
    dplyr::filter(!grepl("INS", var_id)) |>
    separate(var_id, c("chrom", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
    mutate(chrom = as.character(chrom)) |>
    left_join(ld_vcf, join_by(chrom, pos, ref, alt)) |>
    dplyr::select(chrom, pos, rsid, other_allele = ref, effect_allele = alt, r2)

eqtls_ld <-
    eqtls |>
    map(~left_join(., ld_risk_var, join_by(chrom, pos)) |>
            dplyr::select(chrom, pos, rsid, other_allele, effect_allele, p, r2))


loc_list <-
    eqtls_ld |>
    map(~locus(data = .,
               gene = 'IKBKE',
               flank = 5e4,
               LD = "r2",
               ens_db = "EnsDb.Hsapiens.v75"))

# cm_scatter <-
#     gg_scatter(loc_list$cm, 
#                index_snp = 
#                labels = c("rs2297550"),
#                legend_pos = "topright")
# 
# ncm_scatter <-
#     gg_scatter(loc_list$ncm, 
#                labels = "index",
#                legend_pos = "topright")
# 
# nk_scatter <-
#     gg_scatter(loc_list$nk, 
#                labels = "index",
#                legend_pos = "topright")
# 
# pbmc_scatter <-
#     gg_scatter(loc_list$pbmc, 
#                labels = "index",
#                legend_pos = "topright")
# 
# t4_scatter <-
#     gg_scatter(loc_list$t4, 
#                labels = "index",
#                legend_pos = "topright")
# 
# t8_scatter <-
#     gg_scatter(loc_list$t8, 
#                labels = "index",
#                legend_pos = "topright")
# 
# g <- gg_genetracks(loc_list$pbmc)
# 
# p_out <- cm_scatter / ncm_scatter / nk_scatter / pbmc_scatter / t4_scatter / t8_scatter / g
# 
# ggsave("./plots/locuszoom_perez_eqtls.png", p_out, height = 12, width = 6)

gg_scatter_all <- 
    eqtls_ld |>
    bind_rows(.id = "cell") |>
    mutate(r2interval = case_when(r2 >= 0 & r2 < .2 ~ "0.0 - 0.2",
                                  r2 >= .2 & r2 < .4 ~ "0.2 - 0.4",
                                  r2 >= .4 & r2 < .6 ~ "0.4 - 0.6",
                                  r2 >= .6 & r2 < .8 ~ "0.6 - 0.8",
                                  r2 >= .8 & r2 <= 1 ~ "0.8 - 1.0")) |>
    arrange(r2) |>
    mutate(r2interval = fct_rev(r2interval)) |>
    ggplot(aes(x = pos, y = -log10(p))) +
    geom_point(aes(fill = r2interval), shape = 21, stroke = .25, size = 3) +
    scale_x_continuous(limits = c(206600000, 206720000),
                       labels = function(x) x/1e6)+
    scale_fill_manual(values = c("0.0 - 0.2" = "#486CD9",
                                 "0.2 - 0.4" = "#6BEBEC",
                                 "0.4 - 0.6" = "#5DC83B",
                                 "0.6 - 0.8" = "#F3A83B",
                                 "0.8 - 1.0" = "#EB3223")) +
    facet_wrap(~cell, ncol = 1, scale = "free_y") +
    theme_classic() +
    theme(axis.text = element_text(size = 11, family = "Arial"),
          axis.title = element_text(size = 12, family = "Arial"),
          strip.text = element_text(size = 12, family = "Arial"),
          legend.title = element_text(size = 12, family = "Arial"),
          legend.text = element_text(size = 11, family = "Arial")) +
    labs(x = "Chromosome 1 (Mb)",
         fill = "r2 with\nrs2297550:")

g2 <- gg_genetracks(loc_list$pbmc) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          text = element_text(size = 11, family = "Arial"))


ggsave("./plots/locuszoom_perez_eqtls_2.png", 
       gg_scatter_all / g2 + plot_layout(heights = c(1, 1/7)), 
       height = 10, width = 6)



"./data/perez_et_al/science.abf1970_table_s10.txt" |>
    read_tsv(skip = 1) |>
    filter(gene == "IKZF1") |>
    arrange(desc(H4))
