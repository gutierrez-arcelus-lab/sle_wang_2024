library(tidyverse)
library(glue)
library(httr)
library(jsonlite)
library(extrafont)

# Functions
request_datasets_from_api <- function(study_id = "", quant_method = "", sample_group = "",
				      tissue_id = "", study_label = "", tissue_label = "",
				      condition_label = "") {
  
    size <- 1000 #Page size
    start <- 0 #Page start

    parameter_values <- c(study_id,quant_method,sample_group,tissue_id,study_label, 
		       tissue_label,condition_label)
    
    parameter_names <- c('study_id','quant_method','sample_group','tissue_id',
		      'study_label','tissue_label','condition_label')

    while (TRUE) {
	
	URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={size}&start={start}")

	#Adding defined parameters to the request
	for (i in 1:length(parameter_values)) {
	    par <- parameter_values[i]
	    par_name <- parameter_names[i]
	    if (par != "") URL <- glue("{URL}&{par_name}={par}")
	}

	r <- GET(URL, accept_json())
	cont <- content(r, "text", encoding = "UTF-8")

	# If the request was unsuccessful
	if (status_code(r) != 200) {
	    
	    #If we get no results at all, print error
	    if (start == 0) {
		print(glue("Error {status_code(r)}"))
		print(cont)
		return ()
	    }
	  #else just break
	  break
	}

	cont_df <- fromJSON(cont)

	if (start == 0) {
	    responses <- cont_df
	
	} else {
	  responses <- rbind(responses, cont_df)
	}
	
	start <- start + size
    }
    
    return(as_tibble(responses))
}

get_assoc_over_datasets <- function(datasets, variant, gene_id) {
    
    size <- 1000
    first <- TRUE
    
    for (i in rownames(datasets)) {
	row <- datasets[i, ]
	dataset_id <- row$dataset_id

	URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&variant={variant}&gene_id={gene_id}")

	r <- GET(URL, accept_json())

	if (status_code(r) != 200) next

	cont <- content(r, "text", encoding = "UTF-8")
	cont_df <- fromJSON(cont)
	cont_with_metadata <- cbind(cont_df, row)

	if (first) {
	    final_df <- cont_with_metadata
	    first = FALSE
	} else {
	    final_df <- rbind(final_df, cont_with_metadata)
	}
    }
    
    return(as_tibble(final_df))
}


# Data
tissues <- c("macrophage", "monocyte", "neutrophil", "CD4+ T cell", "Treg memory", "LCL", "T cell",
	     "blood", "Tfh cell", "Th17 cell", "Th1 cell", "Th2 cell", "Treg naive", "B cell", 
	     "CD8+ T cell", "CD16+ monocyte", "NK cell")

variant <- "chr1_206470429_C_G"
gene_id <- "ENSG00000263528" 

datasets_ge <- request_datasets_from_api(quant_method = "ge")

datasets_ge_immuno <- filter(datasets_ge, tissue_label %in% tissues) 

associations <- get_assoc_over_datasets(datasets = datasets_ge_immuno, 
					variant = variant, 
					gene_id = gene_id)

plot_df <- 
    associations |>
    select(rsid, beta, study_label, tissue_label, condition_label) |>
    mutate(study_label = sub("_", " (", study_label),
	   study_label = ifelse(grepl("\\(\\d+", study_label), paste0(study_label, ")"), study_label),
	   condition_label = sub("_", " ", condition_label),
	   lab = glue("{study_label} {tissue_label} ({condition_label})"),
	   lab = ifelse(grepl("LCL", lab), sub(" \\(naive\\)", "", lab), lab)) |>
    select(rsid, dataset = lab, beta) |>
    arrange(beta) |>
    mutate(dataset = fct_inorder(dataset))

out_plot <- 
    ggplot(plot_df, aes(x = beta, y = dataset)) +
    geom_point(size = 2) +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-.5, 0, .5)) +
    theme_minimal() +
    theme(
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.x = element_line(color = "grey96"),
	  panel.grid.major.y = element_line(color = "grey96"),
	  plot.title = element_text(family = "Arial", size = 11),
	  axis.text.y = element_text(family = "Arial"),
	  axis.text.x = element_text(family = "Arial"),
	  axis.title.x = element_text(family = "Arial"),
	  plot.title.position = "plot",
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(y = NULL, 
	 title = glue("Beta of rs2297550 ({variant}) on gene expression of\nIKBKE ({gene_id}) in eQTL Catalogue studies"))

ggsave("./plots/rs2297550_betas.png", out_plot, width = 5, height = 4)
