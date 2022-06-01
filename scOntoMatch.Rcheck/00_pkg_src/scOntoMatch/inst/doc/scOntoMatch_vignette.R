## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dpi=300)

## ----install------------------------------------------------------------------

## install from source
# install.packages("devtools")
library(devtools)

devtools::install_github("YY-SONG0718/scOntoMatch")

library(scOntoMatch)
library(ontologyIndex)

## ----load data----------------------------------------------------------------
metadata = '../example_data/metadata.tsv'

anno_col = 'cell_ontology_class'
onto_id_col = 'cell_ontology_id'

obo_file = '../example_data/cl-basic.obo'
propagate_relationships = c('is_a', 'part_of')
ont <- ontologyIndex::get_OBO(obo_file, propagate_relationships = propagate_relationships)

## ----load adata---------------------------------------------------------------

adatas = getAdatas(metadata = metadata, sep = "\t")
adatas

## ----ontoMultiMinimal---------------------------------------------------------
adatas_minimal = scOntoMatch::ontoMultiMinimal(adatas = adatas, ont = ont, anno_col = anno_col, onto_id_col = onto_id_col)

## ----re-annotate--------------------------------------------------------------
adatas$TS_lung$obs[[anno_col]] = as.character(adatas$TS_lung$obs[[anno_col]])

## nk cell can certainly be matched
adatas$TS_lung$obs[which(adatas$TS_lung$obs[[anno_col]] == 'nk cell'), anno_col] = 'natural killer cell'

## there are type 1 and type 2 alveolar fibroblast which both belongs to fibroblast of lung

adatas$TS_lung$obs[which(adatas$TS_lung$obs[[anno_col]] == 'alveolar fibroblast'), anno_col] = 'fibroblast of lung'

## capillary aerocyte is a recently discovered new lung-specific cell type that is good to keep it
## Gillich, A., Zhang, F., Farmer, C.G. et al. Capillary cell-type specialization in the alveolus. Nature 586, 785–789 (2020). https://doi.org/10.1038/s41586-020-2822-7


## ----ontoMultiMinimal_new-----------------------------------------------------
adatas_minimal = scOntoMatch::ontoMultiMinimal(adatas = adatas, ont = ont, anno_col = anno_col, onto_id_col = onto_id_col)

## ----plotOntoTree-------------------------------------------------------------


plotOntoTree(ont = ont, 
                          onts = names(getOntologyId(adatas$TM_lung$obs[['cell_ontology_class']], ont = ont)), 
                          ont_query = names(getOntologyId(adatas$TM_lung$obs[['cell_ontology_class']], ont = ont)),
                          plot_ancestors = TRUE,  roots = 'CL:0000548',
                          fontsize=25)

## ----plotOntoTree_two---------------------------------------------------------


plotOntoTree(ont = ont, 
                          onts = names(getOntologyId(adatas$TS_lung$obs[['cell_ontology_class']], ont = ont)), 
                          ont_query = names(getOntologyId(adatas$TS_lung$obs[['cell_ontology_class']], ont = ont)),
                          plot_ancestors = TRUE,  roots = 'CL:0000548',
                          fontsize=25)

## ----plotOntoTree_minimal-----------------------------------------------------

plotOntoTree(ont = ont, 
                          onts = names(getOntologyId(adatas$TM_lung$obs[['cell_ontology_base']], ont = ont)), 
                          ont_query = names(getOntologyId(adatas$TM_lung$obs[['cell_ontology_base']], ont = ont)),
                          plot_ancestors = TRUE,  roots = 'CL:0000548',
                          fontsize=25)

## ----plotOntoTree_minimal_two-------------------------------------------------

plotOntoTree(ont = ont, 
                          onts = names(getOntologyId(adatas$TS_lung$obs[['cell_ontology_base']], ont = ont)), 
                          ont_query = names(getOntologyId(adatas$TS_lung$obs[['cell_ontology_base']], ont = ont)),
                          plot_ancestors = TRUE,  roots = 'CL:0000548',
                          fontsize=25)

## ----ontoMultiMatch-----------------------------------------------------------

## perform ontoMatch on the original tree

adatas_matched = scOntoMatch::ontoMultiMatch(adatas = adatas_minimal, anno_col = 'cell_ontology_base', onto_id_col = onto_id_col, ont = ont)

## -----------------------------------------------------------------------------
adatas_matched


## ----plotMatchedOntoTree------------------------------------------------------

plts = plotMatchedOntoTree(ont = ont, 
                                 adatas = adatas,
                                 anno_col = 'cell_ontology_mapped', 
                                 onto_id_col = onto_id_col,
                                 roots = 'CL:0000548', fontsize=25)

## -----------------------------------------------------------------------------
plts[[1]]

## ----plotMatchedOntoTree_two--------------------------------------------------
plts[[2]]

## ----getOntologyName----------------------------------------------------------
ont <- ontologyIndex::get_OBO(obo_file, propagate_relationships = c('is_a', 'part_of'))

getOntologyName(adatas[[1]]$obs[[onto_id_col]], ont = ont)


## ----getOntologyId------------------------------------------------------------

getOntologyId(adatas[[2]]$obs[[anno_col]], ont = ont)


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

