# scOntoMatch
scOntoMatch is an R package which unifies ontology annotation of scRNA-seq datasets to make them comparable across studies.


Author&Maintainer: Yuyao Song <ysong@ebi.ac.uk>

Large-scale single cell atlases often have curated annotations using a standard ontology system. 
This package aims to match ontology labels between different datasets to make them comparable across studies. 

## Installation

```
## from source

install.packages("devtools")
devtools::install_github("YY-SONG0718/scOntoMatch")
library(scOntoMatch)
library(ontologyIndex)
```

## Usage
First download the ontology `.obo` file from [the OBO foundry](https://obofoundry.org/). Common ontologies include:

  -  [Cell ontology](https://obofoundry.org/ontology/cl.html)
  -  [Uberon multi-species anatomy ontology](https://obofoundry.org/ontology/uberon.html) 
  -  [Drosophila gross anatomy (FBbt)](https://obofoundry.org/ontology/fbbt.html)
  -  [Zebrafish anatomy and development ontology (ZFA)](https://obofoundry.org/ontology/zfa.html)

Refer to [vignette](https://github.com/YY-SONG0718/scOntoMatch/blob/main/vignettes/scOntoMatch_vignette.Rmd) for detailed usage.

Get input ready

*Note* for anndata .h5ad inputs, goto main branch, for seurat object .rds inputs, goto seurat_obj branch
```
adatas = getAdatas(metadata = metadata, sep = "\t")
ont = ontologyIndex::get_OBO(oboFile, propagate_relationships = c('is_a', 'part_of'), )
```

Trim ontology tree to remove redundant terms

```
adatas_minimal = ontoMultiMinimal(adatas, ont = ont, anno_col = "cell_ontology_type", onto_id_col = "cell_ontology_id")
```

Match ontology cross datasets by direct mapping and mapping descendants to ancestor terms. 

```
adatas_matched = ontoMultiMatch(adatas_minimal, ont = ont, anno_col = "cell_ontology_base")
```


This package also provides convenient plotting functions to help comprehend the hierarchy of cell types in any single cell dataset. 

Plot a ontology tree per dataset
```
plotOntoTree(ont, onts, plot_ancestors=TRUE, ont_query, fontsize = 20)
```
Plot a matched ontology tree for all datasets

```
# use 'animal cells' as root
plts = plotMatchedOntoTree(ont = ont, adatas = adatas, anno_col = "cell_ontology_mapped", roots = 'CL:0000548', fontsize=25)
```                                 

The [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) hosted at EBI provides uniformly analysed and annotated scRNA-Seq data across multiple species. 
Datasets with curated ontology labels are all great inputs to this package. scRNA-seq data stored as h5ad files can be downloaded via the [ftp site](http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/). 
Files that has extension `.project.h5ad` can be pass to `ontoMatch` with `anno_col = 'authors_cell_type_-_ontology_labels'`. 


This package imports functions from [ontologyIndex](https://CRAN.R-project.org/package=ontologyIndex), [ontologyPlot](https://CRAN.R-project.org/package=ontologyPlot) and 
[anndata](https://CRAN.R-project.org/package=anndata)
