
#' read in anndata files as a named list of anndata object
#' @name getAdatas
#' @param metadata a metadaat file indicating name, path to adata file
#' @param sep sep of the metadata file
#' @importFrom anndata read_h5ad
#' @importFrom utils read.table
#' @export
#'
getAdatas <- function(metadata, sep ) {
  metadata <- read.table(metadata, sep = sep, col.names = c("name", "path"))
  adatas <- list()
  n_iter <- nrow(metadata) # Number of iterations of the loop

  # Initializes the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=") # Character used to create the bar
  message("start loading adata objects")
  for(i in 1:n_iter) {

    #---------------------
    # Code to be executed
    #---------------------

    adatas[[metadata[i, "name"]]] <- anndata::read_h5ad(metadata[i, "path"])
    #---------------------

    # Sets the progress bar to the current state
    setTxtProgressBar(pb, i)
  }

  close(pb) # Close the connection

  return(adatas)
}

#' make sure ontology names are all translated to ontology ids
#' while warning, consider manual reannotation
#' @name check_ontology_translate
#' @param adata anndata file
#' @param onts ontology ids from translate
#' @param ont ontologyIndex object
#' @param anno_col annotation column in adata$obs that is translated to onts ids
#' @export
#'


check_ontology_translate <- function(adata, onts, ont, anno_col) {
  if (length(onts) != length(levels(factor(adata$obs[[anno_col]])))) {
    message("warning: some cell type annotations do not have corresponding ontology id, consider manual re-annotate")
    message(paste(levels(factor(adata$obs[[anno_col]]))[!tolower(levels(factor(adata$obs[[anno_col]]))) %in% tolower(ont$name)], collapse = ", ", sep = ", "))
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
#' @export


getOntologyId <- function(onto_name, ont) {
  return(ont$name[names(ont$id[ont$name %in% levels(factor(onto_name))])])
}

#' Get a names list of ontology and id by id
#' @name getOntologyName
#' @param onto_id character vector of ontology ids
#' @param ont ontology object
#' @return a named list mapping ontology id and ontology name
#' @export


getOntologyName <- function(onto_id, ont) {
  return(ont$name[names(ont$id[ont$id %in% levels(factor(onto_id))])])
}


#' Get a names list of ontology and id by id
#' @name fill_query
#' @param all all ontology id to plot tree
#' @param query query ontology id to fill
#' @return a color object to fill query in onto_plot
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
