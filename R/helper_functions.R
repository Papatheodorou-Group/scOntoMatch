
#' read in seurat object .rds files as a named list of seurat object
#' @name getSeuratRds
#' @param metadata a metadata file indicating name, path to 'seurat' rds file
#' @param sep sep of the metadata file
#' @importFrom utils read.table setTxtProgressBar txtProgressBar
#' @return a named list contains data name and the corresponding 'seurat' object
#' @examples
#' \dontrun{
#' getSeuratRds(metadata = 'metadata.tsv', sep = '\t')
#' }
#' @export
#'
getSeuratRds <- function(metadata, sep ) {
  metadata <- read.table(metadata, sep = sep, col.names = c("name", "path"))
  obj_list <- list()
  n_iter <- nrow(metadata) # Number of iterations of the loop

  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=") # Character used to create the bar
  message("start loading seurat rds objects")
  for(i in 1:n_iter) {

    #---------------------
    # Code to be executed
    #---------------------

    obj_list[[metadata[i, "name"]]] <- readRDS(metadata[i, "path"])
    #---------------------

    # Sets the progress bar to the current state
    setTxtProgressBar(pb, i)
  }

  close(pb) # Close the connection

  return(obj_list)
}

#' make sure ontology names are all translated to ontology ids
#' while warning, consider manual reannotation
#' @name check_ontology_translate
#' @param obj seurat rds object
#' @param onts ontology ids from translate
#' @param ont ontologyIndex object
#' @param anno_col annotation column in obj@meta.data that is translated to onts ids
#' @return do not return a value but output messages


check_ontology_translate <- function(obj, onts, ont, anno_col) {
  if (length(onts) != length(levels(factor(obj@meta.data[[anno_col]])))) {
    message("warning: some cell type annotations do not have corresponding ontology id, consider manual re-annotate")
    message(paste(levels(factor(obj@meta.data[[anno_col]]))[!tolower(levels(factor(obj@meta.data[[anno_col]]))) %in% tolower(ont$name)], collapse = ", ", sep = ", "))
  } else {
    (
      message("ontology annotation translate to id successful")

    )
  }
}

#' Get a names list of ontology and id by name
#' @name getOntologyId
#' @param onto_name character vector of ontology names
#' @param ont ontology object
#' @return a named list mapping ontology id and ontology name
#' @examples
#' \dontrun{
#' getOntologyId(onto_name = "epithelial cell of lung", ont = ont)
#' }
#' @export


getOntologyId <- function(ont, onto_name) {
  return(ont$name[names(ont$id[ont$name %in% levels(factor(onto_name))])])
}

#' Get a names list of ontology and id by id
#' @name getOntologyName
#' @param onto_id character vector of ontology ids
#' @param ont ontology object
#' @return a named list mapping ontology id and ontology name
#' @examples
#' \dontrun{
#' getOntologyName(onto_id = "epithelial cell of lung", ont = ont)
#' }
#' @export


getOntologyName <- function(ont, onto_id) {
  return(ont$name[names(ont$id[ont$id %in% levels(factor(onto_id))])])
}


#' Get a names list of ontology and id by id
#' @name fill_query
#' @param all all ontology id to plot tree
#' @param query query ontology id to fill
#' @return a color object to fill query in onto_plot
#' @examples
#' \dontrun{
#' fill_query(all = c("CL0000548", "CL:0000066", "CL:0000082"), query = c("CL:0000082"))
#' }
#' @export

fill_query <- function(all, query) {
  color <- c()
  for (term_now in all) {
    if (term_now %in% query) {
      color <- c(color, "mediumaquamarine")
    } else {
      color <- c(color, "mistyrose")
    }
  }
  return(color)
}
