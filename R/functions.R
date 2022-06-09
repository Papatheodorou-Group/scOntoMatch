#' match ontology terms by direct mapping and mapping descendants to ancestors
#' @name getOntoMapping
#' @param ont the ontology object from get_OBO
#' @param onts1 a character vector of ontology id
#' @param onts2 a character vector of ontology id
#' @return a names list for ontology id mapping looks like ontology_id:ontology_id
#' @importFrom ontologyIndex get_ancestors
#' @export

getOntoMapping <- function(ont, onts1, onts2) {

  ## direct matching
  intersection <- intersect(onts1, onts2)
  mappings <- c()
  mappings[intersection] <- intersection
  message(paste0(" intersection terms: ", intersection))

  onts1 <- onts1[!onts1 %in% intersection]
  onts2 <- onts2[!onts2 %in% intersection]

  onts1 <- minimal_set(ont, onts1)
  onts2 <- minimal_set(ont, onts2)

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
#' @return a names list for ontology id mapping looks like ontology_id:ontology_id
#' @importFrom ontologyIndex get_ancestors
#' @export
#'
#'
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
#' @return an anndata object with .obs[["cell_ontology_base"]]
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


#' get the minimal ontology tree of a list of obj objects by reducing descendant terms to ancestor terms
#' return a named list of obj objects with meta.data[["cell_ontology_base"]] storing the reduced ontology annotation
#' @name ontoMultiMinimal
#' @param obj_list a named list of seurat objects
#' @param ont ontologyIndex object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @return a named list of obj objects with $obs[["cell_ontology_base"]]
#' @importFrom ontologyIndex get_ancestors
#' @export
#'
#'
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
#' Match the ontology annotation of several obj files
#' @name ontoMatch
#' @param adata1 one anndata object
#' @param adata2 the other anndata object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @param ont ontology file loaded via get_OBO
#' @return a list of obj files with annotation ontology mapped to each-other in obs[['cell_type_mapped_ontology']]
#' @importFrom ontologyIndex get_OBO
#' @importFrom anndata read_h5ad
#' @export


ontoMatch <- function(adata1, adata2, anno_col, onto_id_col, ont) {
  message("start matching the ontology annotation")
  ad_one <- adata1
  ad_two <- adata2
  message(paste0("adata1 has cell types: ", paste(levels(factor(ad_one@meta.data[[anno_col]])), collapse = ", ")))
  message(paste0("adata2 has cell types: ", paste(levels(factor(ad_two@meta.data[[anno_col]])), collapse = ", ")))

  if (!is.null(ad_one@meta.data[[onto_id_col]]) & !is.null(ad_two@meta.data[[onto_id_col]])) {
    message("use existing ontology id")
    onts1 <- names(ont$name[names(ont$id[ont$id %in% levels(factor(ad_one@meta.data[[onto_id_col]]))])])
    onts2 <- names(ont$name[names(ont$id[ont$id %in% levels(factor(ad_two@meta.data[[onto_id_col]]))])])
  } else {
    message("translate annotation to ontology id")
    onts1 <- names(ont$name[names(ont$id[tolower(ont$name) %in% tolower(levels(factor(ad_one@meta.data[[anno_col]])))])])
    onts2 <- names(ont$name[names(ont$id[tolower(ont$name) %in% tolower(levels(factor(ad_two@meta.data[[anno_col]])))])])
    if (length(onts1) != length(levels(factor(ad_one@meta.data[[anno_col]])))) {
      message("warning: some cell type annotations do not have corresponding ontology id in obj 1, consider manual re-annotate: ")
      message(paste(levels(factor(ad_one@meta.data[[anno_col]]))[!tolower(levels(factor(ad_one@meta.data[[anno_col]]))) %in% tolower(ont$name)], collapse = ", "))
    }
    if (length(onts2) != length(levels(factor(ad_two@meta.data[[anno_col]])))) {
      message("warning: some cell type annotations do not have corresponding ontology id in obj 2, consider manual re-annotate: ")
      message(paste(levels(factor(ad_two@meta.data[[anno_col]]))[!tolower(levels(factor(ad_two@meta.data[[anno_col]]))) %in% tolower(ont$name)], collapse = ", "))
    }
  }


  mappings <- getOntoMapping(ont = ont, onts1 = onts1, onts2 = onts2)


  ad_one@meta.data[["cell_type_mapped_ontology"]] <- as.character(ad_one@meta.data[[anno_col]])
  ad_two@meta.data[["cell_type_mapped_ontology"]] <- as.character(ad_two@meta.data[[anno_col]])

  for (fromTerm in names(mappings)) {
    toTerm <- mappings[fromTerm]
    fromName <- ont$name[names(ont$id[ont$id == fromTerm])]
    toName <- ont$name[names(ont$id[ont$id == toTerm])]
    message(paste("mapping from name: ", fromName, " to name: ", toName, sep = ""))
    ad_one@meta.data[which(tolower(ad_one@meta.data[[anno_col]]) == tolower(fromName)), "cell_type_mapped_ontology"] <- toName
    ad_two@meta.data[which(tolower(ad_two@meta.data[[anno_col]]) == tolower(fromName)), "cell_type_mapped_ontology"] <- toName
  }

  message(paste0("after mapping, adata1 has cell types: ", paste(levels(factor(ad_one@meta.data[["cell_type_mapped_ontology"]])), collapse = ", ")))
  message(paste0("after mapping, adata2 has cell types: ", paste(levels(factor(ad_two@meta.data[["cell_type_mapped_ontology"]])), collapse = ", ")))

  return(list(ad_one, ad_two))
}

#' Core function of scOntoMatch
#' Match the ontology annotation of multiple obj objects
#' @name ontoMultiMatch
#' @param obj_list a namesd list of seurat objects to match
#' @param ont ontologyIndex object
#' @param anno_col the cell ontology text annotation column name
#' @param onto_id_col if also have ontology id column for direct mapping
#' @return a list of obj objects with annotation ontology mapped to each-other in obs[['cell_ontology_mapped']]
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


#' In a matched ontology tree there could be simultaneously ancestral terms and descendant terms within each dataset,
#' this function takes the output of ontoMatch and merge descendant terms to ancestor terms, while not over-merge
#' @name ontoMatchMinimal
#' @param ont ontology object
#' @param obj_list a list of obj files from the output of ontoMatch
#' @return a list of obj files with minimal ontology mapped to each other in obs[['cell_type_mapped_ontology_base']]
#' @importFrom ontologyIndex get_OBO get_ancestors get_descendants minimal_set
#' @export


ontoMatchMinimal <- function(obj_list, ont) {

  # the new ontology terms that has been mapped between datasets
  new_all <- unique(c(names(getOntologyId(obj_list[[1]]@meta.data[["cell_type_mapped_ontology"]], ont = ont)), names(getOntologyId(obj_list[[2]]@meta.data[["cell_type_mapped_ontology"]], ont = ont))))

  # minimal_set(ont, new_all) is all the leaf terms in the current joint ontology tree
  # setdiff(new_all, minimal_set(ont, new_all)) contains intermediate terms
  # minimal_set(ont, setdiff(new_all, minimal_set(ont, new_all))) contains leaf terms among intermediate terms to avoid over-merging
  # that gives the removed_terms for us to match to its ancestor in intermediate terms
  removed_terms_all <- setdiff(new_all, minimal_set(ont, new_all))
  removed_terms <- setdiff(removed_terms_all, minimal_set(ont, removed_terms_all))

  new_terms_one <- names(getOntologyId(obj_list[[1]]@meta.data[["cell_type_mapped_ontology"]], ont = ont))
  new_terms_two <- names(getOntologyId(obj_list[[2]]@meta.data[["cell_type_mapped_ontology"]], ont = ont))

  # get the mapping between elements to match to ancestor and the ancestor
  # it gives the same result to use new_terms_one or new_terms_two in the following line
  to_minimize <- new_terms_one[new_terms_one %in% get_descendants(ontology = ont, roots = removed_terms) & !(new_terms_one %in% removed_terms)]
  common_base <- lapply(structure(to_minimize, names = to_minimize), function(x) removed_terms[unlist(lapply(removed_terms, function(y) x %in% get_descendants(ontology = ont, roots = y)))])
  common_base <- common_base[unlist(lapply(common_base, function(x) length(x) > 0))]
  mappings <- unlist(common_base, use.names = TRUE)

  # from previously unified annotation
  obj_list[[1]]@meta.data[["cell_type_mapped_ontology_base"]] <- obj_list[[1]]@meta.data[["cell_type_mapped_ontology"]]
  obj_list[[2]]@meta.data[["cell_type_mapped_ontology_base"]] <- obj_list[[2]]@meta.data[["cell_type_mapped_ontology"]]

  ## perform mapping
  for (fromTerm in names(mappings)) {
    toTerm <- mappings[fromTerm]
    fromName <- ont$name[names(ont$id[ont$id == fromTerm])]
    toName <- ont$name[names(ont$id[ont$id == toTerm])]
    message(paste("mapping from name: ", fromName, " to name: ", toName, sep = ""))
    obj_list[[1]]@meta.data[which(obj_list[[1]]@meta.data[["cell_type_mapped_ontology"]] == fromName), "cell_type_mapped_ontology_base"] <- toName
    obj_list[[2]]@meta.data[which(obj_list[[2]]@meta.data[["cell_type_mapped_ontology"]] == fromName), "cell_type_mapped_ontology_base"] <- toName
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
#' @importFrom ontologyPlot onto_plot
#' @importFrom purrr flatten_chr
#' @importFrom ontologyIndex get_ancestors intersection_with_descendants
#' @export
#'
plotMatchedOntoTree <- function(obj_list, ont, anno_col, onto_id_col, roots = c("CL:0000548"), ...) {
  plots <- list()
  onts <- suppressMessages(ontoTranslate(obj_list, ont, onto_id_col, anno_col = anno_col))
  all <- unique(flatten_chr(onts))

  for (i in seq(1, length(obj_list))) {
    plots[[names(obj_list[i])]] <- plotOntoTree(ont = ont, onts = all, ont_query = names(getOntologyId(obj_list[[i]]@meta.data[[anno_col]], ont = ont)), plot_ancestors = TRUE, roots = roots, ...)
  }

  return(plots)
}
