library(tidyverse)
library(httr)
library(glue)
library(jsonlite)

# Function to request data in a region
request_associations <- function(dataset_id, chromosome_id, position, distance) {
    
    page_size <- 1000
    page_start <- 0
    range_start <- position - distance
    range_end <- position + distance

    while (TRUE) {
	
	URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={page_size}&start={page_start}&pos={chromosome_id}:{range_start}-{range_end}")
	r <- GET(URL, accept_json())

	# If Error == 500, try again up to 30 times
	ctn <- 0
	while (status_code(r) == 500 && ctn < 30) {
	    r <- GET(URL, accept_json())
	    ctn <- ctn + 1L
	    Sys.sleep(5)
	}

	cont <- content(r, "text", encoding = "UTF-8")

	# If the request was unsuccessful
	if (status_code(r) != 200) {
	    
	    #If we get no results at all, print error
	    if (page_start == 0) {
		stop(glue("{status_code(r)}\n{cont}"))
	    }
	    
	    #else just break
	    break
	}

	cont_df <- fromJSON(cont)
	
	if (page_start == 0) {
	    responses <- cont_df
	} else {
	    responses <- bind_rows(responses, cont_df)
	}

	page_start <- page_start + page_size
    }
    
    return(responses)
}

safe_request <- safely(.f = request_associations)

# Region of interest
chr <- 7
pos <- 50266267 
window_size <- 5e5 
out <- "/lab-share/IM-Gutierrez-e2/Public/vitor/colabs/sle_ikke/ikzf1/results"

# eQTL catalogue dataset
query_id <- commandArgs(TRUE)[1]

# Request eQTL catalogue data
associations <- 
    safe_request(dataset_id = query_id,
		 chromosome_id = chr,
		 position = pos,
		 distance = window_size)

if (is.data.frame(associations$result)) {
    
    out_file <-
	file.path(out, "%s.tsv") |>
	sprintf(query_id)

    write_tsv(associations$result, out_file)
}

if (!is.null(associations$error)) {

    out_err <-
	file.path(out, "%s.err") |>
	sprintf(query_id)

    write_lines(associations$error, out_err)
}
