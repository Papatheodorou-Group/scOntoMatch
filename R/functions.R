#' match ontology terms by direct mapping and mapping descendants to ancestors
#' @name getOntoMapping
#' @param ont the ontology object from get_OBO
#' @param onts1 a character vector of ontology id
#' @param onts2 a character vector of ontology id
#' @return a names list for ontology id mapping looks like ontology_id:ontology_id
#' @importFrom ontologyIndex get_ancestors
#' @export

getOntoMapping <- function(ont, onts1, onts2){

  intersection <- intersect(onts1, onts2)
  mappings = c()
  mappings[intersection] = intersection
  message(paste0("intersection terms: ", intersection))

  onts1 <- onts1[! onts1 %in% intersection]
  onts2 <- onts2[! onts2 %in% intersection]

  ancestors_onts_1 <- lapply(onts1, function(x) ontologyIndex::get_ancestors(ont, x))
  ancestors_onts_2 <- lapply(onts2, function(x) ontologyIndex::get_ancestors(ont, x))

  one_ancestors_of_two <- lapply(structure(onts1, names=onts1), function(x) onts2[unlist(lapply(ancestors_onts_2, function(y) x %in% y))])
  one_ancestors_of_two <- one_ancestors_of_two[unlist(lapply(one_ancestors_of_two, function(x) length(x) > 0))]

  for (ancestor in names(one_ancestors_of_two)){
    mappings[one_ancestors_of_two[[ancestor]]] <- ancestor
    onts1 <- onts1[onts1 != ancestor]
    onts2 <- onts2[! onts2 %in% one_ancestors_of_two[[ancestor]]]
  }


  ancestors_onts_1 <- lapply(onts1, function(x) ontologyIndex::get_ancestors(ont, x))
  ancestors_onts_2 <- lapply(onts2, function(x) ontologyIndex::get_ancestors(ont, x))

  two_ancestors_of_one <- lapply(structure(onts2, names = onts2), function(x) onts1[unlist(lapply(ancestors_onts_1, function(y) x %in% y))])
  two_ancestors_of_one <- two_ancestors_of_one[unlist(lapply(two_ancestors_of_one, function(x) length(x) > 0))]

  for (ancestor in names(two_ancestors_of_one)){
    mappings[two_ancestors_of_one[[ancestor]]] <- ancestor
    onts2 <- onts2[onts2 != ancestor]
    onts1 <- onts1[! onts1 %in% two_ancestors_of_one[[ancestor]]]
  }
  return(mappings)

}



#' Core function of scOntoMatch
#' Match the ontology annotation of several adata files
#' @name ontoMatch
#' @param adata1 one anndata file path
#' @param adata2 the other anndata file path
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @param obo_file obo file path
#' @param propagate_relationships relationships for reading in obo file
#' @return a list of adata files with annotation ontology mapped to each-other
#' @importFrom ontologyIndex get_OBO
#' @importFrom anndata read_h5ad write_h5ad
#' @export


ontoMatch <- function(adata1, adata2, anno_col='authors_cell_type_-_ontology_labels', onto_id_col='authors_cell_type_-_ontology_labels_ontology', obo_file, propagate_relationships = c('is_a', 'part_of'), ...) {
  message("start matching the ontology annotation")
  ont <- ontologyIndex::get_OBO(obo_file, propagate_relationships = propagate_relationships, ...)
  ad_one = anndata::read_h5ad(adata1)
  ad_two = anndata::read_h5ad(adata2)
  message(paste0("adata1 has cell types: ", paste(levels(factor(ad_one$obs[[anno_col]])), collapse = ', ')))
  message(paste0("adata2 has cell types: ", paste(levels(factor(ad_two$obs[[anno_col]])), collapse = ', ')))

  if(!is.null(ad_one$obs[[onto_id_col]]) & !is.null(ad_two$obs[[onto_id_col]])) {
    message("use existing ontology id")
    onts1=names(ont$name[names(ont$id[ont$id %in% levels(factor(ad_one$obs[[onto_id_col]]))])])
    onts2=names(ont$name[names(ont$id[ont$id %in% levels(factor(ad_two$obs[[onto_id_col]]))])])

  } else {
    message("translate annotation to ontology id")
    onts1=names(ont$name[names(ont$id[ont$name %in% levels(factor(ad_one$obs[[anno_col]]))])])
    onts2=names(ont$name[names(ont$id[ont$name %in% levels(factor(ad_two$obs[[anno_col]]))])])
    if(length(onts1) != length(levels(factor(ad_one$obs[[onto_id_col]])))) {
      message("warning: some cell type annotations do not have corresponding ontology id in adata 1, consider manual re-annotate")
      message(paste(levels(factor(ad_one$obs[[anno_col]]))[!levels(factor(ad_one$obs[[anno_col]])) %in% ont$name], collapse = ', '))
    }
    if(length(onts2) != length(levels(factor(ad_two$obs[[onto_id_col]])))) {
      message("warning: some cell type annotations do not have corresponding ontology id in adata 2, consider manual re-annotate")
      message(paste(levels(factor(ad_two$obs[[anno_col]]))[!levels(factor(ad_two$obs[[anno_col]])) %in% ont$name], collapse = ', '))
    }
  }

  mappings = getOntoMapping(ont = ont, onts1 = onts1, onts2 = onts2)

  ad_one$obs[["cell_type_mapped_ontology"]] = as.character(ad_one$obs[[anno_col]])
  ad_two$obs[["cell_type_mapped_ontology"]] = as.character(ad_two$obs[[anno_col]])

  for (fromTerm in names(mappings)){
    toTerm <- mappings[fromTerm]
    fromName = ont$name[names(ont$id[ont$id == fromTerm])]
    toName = ont$name[names(ont$id[ont$id == toTerm])]
    ad_one$obs[which(ad_one$obs[[anno_col]] == fromName), "cell_type_mapped_ontology"] <- toName
    ad_two$obs[which(ad_two$obs[[anno_col]] == fromName), "cell_type_mapped_ontology"] <- toName

  }

  message(paste0("after mapping, adata1 has cell types: ", paste(levels(factor(ad_one$obs[["cell_type_mapped_ontology"]])), collapse = ', ')))
  message(paste0("after mapping, adata2 has cell types: ", paste(levels(factor(ad_two$obs[["cell_type_mapped_ontology"]])), collapse = ', ')))

  return(list(ad_one, ad_two))


}


#' Get a names list of ontology and id by name
#' @name getOntologyId
#' @param adata adata object
#' @param anno_col the cell ontology text annotation column name
#' @param obo_object obo object
#' @return a named list mapping ontology id and ontology name
#' @export


getOntologyId <- function(adata, anno_col='authors_cell_type_-_ontology_labels', obo_object){

  return(obo_object$name[names(obo_object$id[obo_object$name %in% levels(factor(adata$obs[[anno_col]]))])])

}

#' Get a names list of ontology and id by id
#' @name getOntologyName
#' @param adata adata object
#' @param onto_id_col if also have ontology id column for direct mapping
#' @param obo_object obo object
#' @return a named list mapping ontology id and ontology name
#' @export


getOntologyName <- function(adata, onto_id_col='authors_cell_type_-_ontology_labels_ontology', obo_object){

    return(obo_object$name[names(obo_object$id[obo_object$id %in% levels(factor(adata$obs[[onto_id_col]]))])])
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


#' Plot a tree representation of ontology terms
#' @name plotOntoTree
#' @param ont ontology object
#' @param onts ontology ids to plot
#' @param plot_ancestors if plot ancestors or not
#' @param ont_query query ontology to highlight in the tree
#' @return an ontology tree plot
#' @importFrom ontologyPlot onto_plot
#' @importFrom ontologyIndex get_ancestors
#' @export

plotOntoTree <- function(ont, onts, plot_ancestors=TRUE, ont_query=NULL, ...){

  if(plot_ancestors){

    terms = ontologyIndex::get_ancestors(ont, onts)
  } else {

    terms = onts
  }

  plt = ontologyPlot::onto_plot(ontology = ont,
                                terms=terms,
                                fillcolor = fill_query(all = terms, query = ont_query), ...)
  return(plt)

}
