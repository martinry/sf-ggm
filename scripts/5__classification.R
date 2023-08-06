library(httr)
library(jsonlite)

panther.classify <- function(accessions) {
    annotations <- vector("list", length(accessions))
    
    for (i in seq_along(accessions)) {
        uniprot_accession <- accessions[i]
        organism <- "9606"
        print(paste0("Processing index: ", i, ", accession: ", uniprot_accession))
        
        url <- paste0("https://pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList=", uniprot_accession, "&organism=", organism)
        response <- GET(url)
        
        
        
        
        # Check if the request was successful
        if (response$status_code == 200) {
            content <- httr::content(response, "text")
            
            # Convert the JSON text to a list
            data <- fromJSON(content)
            
            if("mapped_genes" %in% names(data$search)) {
                
                if(!is.null(data$search$mapped_genes$gene$annotation_type_list$annotation_data_type$content)) {
                    
                    # Which column holds protein classification?
                    pc = data$search$mapped_genes$gene$annotation_type_list$annotation_data_type$content %>% grep("ANNOT_TYPE_ID_PANTHER_PC", .)
                    
                }
                    
                    anno = tryCatch(
                        expr = {
                            data$search$mapped_genes$gene$annotation_type_list$annotation_data_type$annotation_list$annotation[[pc]]$name
                        },
                        error = function(e){
                            print(e)
                            e
                        }
                    )
                
                if(inherits(anno, "error")) {
                    annotations[[i]] = NA
                } else {
                    annotations[[i]] = anno
                }
            } else {
                annotations[[i]] = NA
            }
            
        } else {
            print(paste0("Failed to fetch data from PantherDB for accession: ", uniprot_accession))
            annotations[[i]] <- NA  # Store NA if there is a failure
        }
    }
    
    return(unlist(annotations))
}
