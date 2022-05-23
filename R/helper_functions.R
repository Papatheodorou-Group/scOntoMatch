
#' read in anndata files as a named list of anndata object
#' @name getAdatas
#' @param metadata a metadaat file indicating name, path to adata file
#' @param sep sep of the metadata file, default "\t"
#' @importFrom anndata read_h5ad
#' @export
#'
getAdatas = function(metadata, sep = "\t"){

  metadata = read.table(metadata, sep = sep, col.names = c("name", "path"))
  adatas = list()
  for(i in seq(1, nrow(metadata))){

    adatas[[metadata[i, 'name']]] = anndata::read_h5ad(metadata[i, 'path'])
  }

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


check_ontology_translate = function(adata, onts, ont, anno_col){

  if(length(onts) != length(levels(factor(adata$obs[[anno_col]])))) {
    message("warning: some cell type annotations do not have corresponding ontology id, consider manual re-annotate")
    message(paste(as.character(levels(factor(adata$obs[[anno_col]]))[!tolower(levels(factor(adata$obs[[anno_col]]))) %in% tolower(ont$name)], collapse = ', ')))
  } else (
    message("ontology annotation translate to id successful")

  )

}

#' Get a names list of ontology and id by name
#' @name getOntologyId
#' @param onto_name character vector of ontology names
#' @param ont ontology object
#' @return a named list mapping ontology id and ontology name
#' @export


getOntologyId <- function(onto_name, ont){

  return(ont$name[names(ont$id[ont$name %in% levels(factor(onto_name))])])

}

#' Get a names list of ontology and id by id
#' @name getOntologyName
#' @param onto_id character vector of ontology ids
#' @param ont ontology object
#' @return a named list mapping ontology id and ontology name
#' @export


getOntologyName <- function(onto_id, ont){

  return(ont$name[names(ont$id[ont$id %in% levels(factor(onto_id))])])
}


#' Helper function to fill queries in plotOntoTree
fill_query = function(all, query) {

  color = c()
  for(term_now in all){

    if(term_now %in% query) {
      color = c(color, "mediumaquamarine")
    } else {
      color = c(color, 'mistyrose')
    }
  }
  return(color)
}