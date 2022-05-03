#'
#' @export
setGeneric(
  name = "ontoMatch",
  def = function(adatas,
                 anno_col,
                 obo_file,
                 propagate_relationships = c('is_a', 'part_of')){
    standardGeneric("ontoMatch")
  }
)

