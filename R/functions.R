#' match ontology terms by direct mapping and mapping descendants to ancestors
#' @name getOntoMapping
#' @param ont the ontology object from get_OBO
#' @param onts1 a character vector of ontology id
#' @param onts2 a character vector of ontology id
#' @return a named list for ontology id mapping looks like ontology_id:ontology_id
#' @examples
#' \dontrun{
#' getOntoMapping(ont = ont, onts1 = "CL:0000548", onts2 = c("CL0000548", "CL:0000066"))
#' }
#' @importFrom ontologyIndex get_ancestors minimal_set
#' @export

getOntoMapping <- function(ont, onts1, onts2) {

  ## direct matching
  intersection <- intersect(onts1, onts2)
  mappings <- c()
  mappings[intersection] <- intersection
  message(paste0(" intersection terms: ", intersection))

  onts1 <- onts1[!onts1 %in% intersection]
  onts2 <- onts2[!onts2 %in% intersection]

  onts1 <- ontologyIndex::minimal_set(ont, onts1)
  onts2 <- ontologyIndex::minimal_set(ont, onts2)

  ancestors_onts_1 <- lapply(onts1, function(x) ontologyIndex::get_ancestors(ont, x))
  ancestors_onts_2 <- lapply(onts2, function(x) ontologyIndex::get_ancestors(ont, x))

  # a term in onts1 is the ancestor of other term(s) in onts2
  one_ancestors_of_two <- lapply(structure(onts1, names = onts1), function(x) onts2[unlist(lapply(ancestors_onts_2, function(y) x %in% y))])
  one_ancestors_of_two <- one_ancestors_of_two[unlist(lapply(one_ancestors_of_two, function(x) length(x) > 0))]

  for (ancestor in names(one_ancestors_of_two)) {
    mappings[one_ancestors_of_two[[ancestor]]] <- ancestor
    onts1 <- onts1[onts1 != ancestor]
    onts2 <- onts2[!onts2 %in% one_ancestors_of_two[[ancestor]]]
  }


  ancestors_onts_1 <- lapply(onts1, function(x) ontologyIndex::get_ancestors(ont, x))
  ancestors_onts_2 <- lapply(onts2, function(x) ontologyIndex::get_ancestors(ont, x))

  # a term in onts2 is the ancestor of other term(s) in onts1

  two_ancestors_of_one <- lapply(structure(onts2, names = onts2), function(x) onts1[unlist(lapply(ancestors_onts_1, function(y) x %in% y))])
  two_ancestors_of_one <- two_ancestors_of_one[unlist(lapply(two_ancestors_of_one, function(x) length(x) > 0))]

  for (ancestor in names(two_ancestors_of_one)) {
    mappings[two_ancestors_of_one[[ancestor]]] <- ancestor
    onts2 <- onts2[onts2 != ancestor]
    onts1 <- onts1[!onts1 %in% two_ancestors_of_one[[ancestor]]]
  }

  return(mappings)
}


#' match descendant terms to ancestor terms within a dataset
#' @name getOntoMinimal
#' @param ont the ontology object from get_OBO
#' @param onts a character vector of ontology id
#' @return a named list for ontology id mapping looks like ontology_id:ontology_id
#' @examples
#' \dontrun{
#' getOntoMinimal(ont = ont, onts = c("CL0000548", "CL:0000066", "CL:0000082"))
#' }
#' @importFrom ontologyIndex get_ancestors
#' @export
#'

getOntoMinimal <- function(ont, onts) {
  ansc_to_desc <- lapply(structure(onts, names = onts), function(x) onts[unlist(lapply(onts, function(y) x %in% get_ancestors(ont, y) & x != y))])
  ansc_to_desc <- ansc_to_desc[unlist(lapply(ansc_to_desc, function(x) length(x) > 0))]

  desc_to_ansc <- list()

  for (ancestor in names(ansc_to_desc)) {
    desc_to_ansc[ansc_to_desc[[ancestor]]] <- ancestor
  }

  return(desc_to_ansc)
}


#' get the minimal ontology tree of a dataset by reducing descendant terms to ancestor terms
#' return  obj meta.data[["cell_ontology_base"]] storing the reduced ontology annotation
#' @name ontoMinimal
#' @param obj the seurat object
#' @param ont ontologyIndex object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @return an seurat object with  meta.data[["cell_ontology_base"]]
#' @examples
#' \dontrun{
#' ontoMinimal(obj = seurat_obj, ont = ont, anno_col = 'ontology_name', onto_id_col = 'ontology_id')
#' }
#' @importFrom ontologyIndex get_ancestors
#' @export
#'

ontoMinimal <- function(obj, ont, anno_col, onto_id_col) {
  if (!is.null(obj@meta.data[[onto_id_col]])) {
    message("use existing ontology id")
    onts1 <- names(ont$name[names(ont$id[ont$id %in% levels(factor(obj@meta.data[[onto_id_col]]))])])
  } else {
    message("translate annotation to ontology id")
    onts1 <- names(ont$name[names(ont$id[tolower(ont$name) %in% tolower(levels(factor(obj@meta.data[[anno_col]])))])])

    if (length(onts1) != length(levels(factor(obj@meta.data[[onto_id_col]])))) {
      message("warning: some cell type annotations do not have corresponding ontology id in obj, consider manual re-annotate: ")
      message(paste(levels(factor(obj@meta.data[[anno_col]]))[!tolower(levels(factor(obj@meta.data[[anno_col]]))) %in% tolower(ont$name)], collapse = ", "))
    }
  }

  desc_to_ansc <- getOntoMinimal(ont = ont, onts = onts1)

  obj@meta.data[["cell_ontology_base"]] <- as.character(obj@meta.data[[anno_col]])

  for (fromTerm in names(desc_to_ansc)) {
    toTerm <- desc_to_ansc[fromTerm]
    fromName <- ont$name[names(ont$id[ont$id == fromTerm])]
    toName <- ont$name[names(ont$id[ont$id == toTerm])]
    message(paste("mapping from name: ", fromName, " to name: ", toName, sep = ""))
    obj@meta.data[which(tolower(obj@meta.data[[anno_col]]) == tolower(fromName)), "cell_ontology_base"] <- toName
  }
  message(paste0("after matching to base level ontology, obj has cell types: ", paste(levels(factor(obj@meta.data[["cell_ontology_base"]])), collapse = ", ")))
  return(obj)
}



#' translate named list of obj_list to named list of cell ontology ids per obj
#' @name ontoTranslate
#' @param obj_list a named list of seurat object
#' @param ont ontologyIndex object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @return a named list of cell ontology ids
#' @examples
#' \dontrun{
#' ontoTranslate(seurat_obj_list, ont, 'ontology_name', 'ontology_id')
#' }
#' @export
#'

ontoTranslate <- function(obj_list, ont, onto_id_col, anno_col) {
  onts <- list()

  if (all(unlist(lapply(obj_list, function(x) !is.null(x@meta.data[[onto_id_col]]))))) {
    message("use existing ontology id")
    for (i in seq(1, length(obj_list))) {
      onts[[names(obj_list[i])]] <- names(ont$name[names(ont$id[ont$id %in% levels(factor(obj_list[[i]]@meta.data[[onto_id_col]]))])])
    }
  } else {
    message("translate annotation to ontology id")

    for (i in seq(1, length(obj_list))) {
      message(paste0("translating ", names(obj_list[i])))
      onts[[names(obj_list[i])]] <- names(ont$name[names(ont$id[tolower(ont$name) %in% tolower(levels(factor(obj_list[[i]]@meta.data[[anno_col]])))])])
      check_ontology_translate(obj = obj_list[[i]], onts = onts[[names(obj_list[i])]], ont = ont, anno_col = anno_col)
    }
  }

  return(onts)
}


#' get the minimal ontology tree of a list of seurat objects by reducing descendant terms to ancestor terms
#' return a named list of seurat objects with meta.data[["cell_ontology_base"]] storing the reduced ontology annotation
#' @name ontoMultiMinimal
#' @param obj_list a named list of seurat objects
#' @param ont ontologyIndex object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @return a named list of seurat objects with meta.data[["cell_ontology_base"]]
#' @examples
#' \dontrun{
#' ontoMultiMinimal(seurat_obj_list, ont, "cell_ontology_base", 'ontology_id')
#' }
#' @importFrom ontologyIndex get_ancestors
#' @export
#'

ontoMultiMinimal <- function(obj_list, ont, anno_col = "cell_ontology_base", onto_id_col) {
  onts <- ontoTranslate(obj_list = obj_list, ont = ont, onto_id_col = onto_id_col, anno_col = anno_col)

  for (i in seq(1, length(obj_list))) {
    desc_to_ansc <- getOntoMinimal(ont = ont, onts = onts[[names(obj_list[i])]])
    obj_list[[i]]@meta.data[["cell_ontology_base"]] <- as.character(obj_list[[i]]@meta.data[[anno_col]])
    for (fromTerm in names(desc_to_ansc)) {
      toTerm <- desc_to_ansc[fromTerm]
      fromName <- ont$name[names(ont$id[ont$id == fromTerm])]
      toName <- ont$name[names(ont$id[ont$id == toTerm])]
      message(paste("mapping from name: ", fromName, " to name: ", toName, sep = ""))
      obj_list[[i]]@meta.data[which(tolower(obj_list[[i]]@meta.data[[anno_col]]) == tolower(fromName)), "cell_ontology_base"] <- toName
    }
    message(paste0("after matching to base level ontology, ", names(obj_list[i]), " has cell types: "))
    message(paste(as.character(levels(factor(obj_list[[i]]@meta.data[["cell_ontology_base"]]))), sep = ", ", collapse = ", "))
  }

  return(obj_list)
}



#' Match descendants to ancestors in multiple ontology id lists
#' @name getOntoMultiMapping
#' @param onts named list of ontology ids
#' @param ont ontologyIndex object
#' @return a named character of mapping from:mapping to
#' @examples
#' \dontrun{
#' getOntoMultiMapping(ont = ont, onts = c("CL0000548", "CL:0000066", "CL:0000082"))
#' }
#' @importFrom purrr flatten_chr
#' @export

getOntoMultiMapping <- function(ont, onts) {

  ## direct matching
  intersection <- Reduce(intersect, onts)
  mappings <- c()
  mappings[intersection] <- intersection
  message(paste0("intersection terms: ", intersection))

  onts_new <- lapply(onts, FUN = function(x) x[!x %in% intersection])
  onts_all <- flatten_chr(onts_new)

  desc_to_ansc <- getOntoMinimal(ont = ont, onts = onts_all)
  if(length(desc_to_ansc) == 0){
    return(mappings)

  } else {
    for (i in seq(1, length(desc_to_ansc))) {
      mappings[names(desc_to_ansc[i])] <- desc_to_ansc[[i]]
    }
  }

  return(mappings)
}



#' Core function of scOntoMatch
#' Match the ontology annotation of multiple seurat objects
#' @name ontoMultiMatch
#' @param obj_list a namesd list of seurat objects to match
#' @param ont ontologyIndex object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @return a list of seurat objects with annotation ontology mapped to each-other in obs[['cell_ontology_mapped']]
#' @examples
#' \dontrun{
#' ontoMultiMatch(seurat_obj_list, ont, "ontology_name", 'ontology_id')
#' }
#' @export

ontoMultiMatch <- function(obj_list, anno_col, onto_id_col, ont) {
  onts <- ontoTranslate(obj_list = obj_list, ont = ont, onto_id_col = onto_id_col, anno_col = anno_col)
  mappings <- getOntoMultiMapping(ont, onts = onts)

  for (i in seq(1, length(obj_list))) {
    message(paste0("processing ", names(obj_list[i])))
    obj_list[[i]]@meta.data[["cell_ontology_mapped"]] <- as.character(obj_list[[i]]@meta.data[[anno_col]])
    for (fromTerm in names(mappings)) {
      toTerm <- mappings[fromTerm]
      fromName <- ont$name[names(ont$id[ont$id == fromTerm])]
      toName <- ont$name[names(ont$id[ont$id == toTerm])]
      message(paste("mapping from name: ", fromName, " to name: ", toName, sep = ""))
      obj_list[[i]]@meta.data[which(tolower(obj_list[[i]]@meta.data[[anno_col]]) == tolower(fromName)), "cell_ontology_mapped"] <- toName
    }
    message(paste0("after matching across datasets, ", names(obj_list[i]), " has cell types: "))
    message(paste(as.character(levels(factor(obj_list[[i]]@meta.data[["cell_ontology_mapped"]]))), sep = ", ", collapse = ", "))
  }

  return(obj_list)
}




#' Plot a tree representation of ontology terms
#' @name plotOntoTree
#' @param ont ontology object
#' @param onts ontology ids to plot
#' @param plot_ancestors if plot ancestors or not
#' @param ont_query query ontology to highlight in the tree
#' @param ... additional parameters for ontologyPlot::onto_plot
#' @param roots root ontology in tree, default "animal cells" in cell ontology
#' @return an ontology tree plot
#' @examples
#' \dontrun{
#' plotOntoTree(ont = ont, onts = c("CL:0000066", "CL:0000082"), ont_query = c("CL:0000082"))
#' }
#' @importFrom ontologyPlot onto_plot
#' @importFrom ontologyIndex get_ancestors intersection_with_descendants
#' @export

plotOntoTree <- function(ont, onts, plot_ancestors = TRUE, ont_query = NULL, roots = c("CL:0000548"), ...) {
  if (plot_ancestors) {
    terms <- ontologyIndex::get_ancestors(ont, onts)
  } else {
    terms <- onts
  }

  if (!is.null(roots)) {
    terms <- ontologyIndex::intersection_with_descendants(ontology = ont, roots = roots, terms = terms)
  }

  plt <- ontologyPlot::onto_plot(
    ontology = ont,
    terms = terms,
    fillcolor = fill_query(all = terms, query = ont_query), ...
  )
  return(plt)
}


#' Plot a ontology tree with matched ontology from ontoMatch
#' @name plotMatchedOntoTree
#' @param ont ontology object
#' @param obj_list a list of seurat obj files as the output of ontoMatch
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @param ... additional parameters for ontologyPlot::onto_plot
#' @param roots root ontology in tree to plot, default "animal cells" in cell ontology
#' @return a lit of matched ontology tree plot
#' @examples
#' \dontrun{
#' plotMatchedOntoTree(seurat_obj_list, ont, 'cell_ontology_mapped', 'ontology_id')
#' }
#' @importFrom ontologyPlot onto_plot
#' @importFrom purrr flatten_chr
#' @importFrom ontologyIndex get_ancestors intersection_with_descendants
#' @export
#'
plotMatchedOntoTree <- function(obj_list, ont, anno_col = 'cell_ontology_mapped', onto_id_col, roots = c("CL:0000548"), ...) {
  plots <- list()
  onts <- suppressMessages(ontoTranslate(obj_list, ont, onto_id_col, anno_col = anno_col))
  all <- unique(flatten_chr(onts))

  for (i in seq(1, length(obj_list))) {
    plots[[names(obj_list[i])]] <- plotOntoTree(ont = ont, onts = all, ont_query = names(getOntologyId(obj_list[[i]]@meta.data[[anno_col]], ont = ont)), plot_ancestors = TRUE, roots = roots, ...)
  }

  return(plots)
}
