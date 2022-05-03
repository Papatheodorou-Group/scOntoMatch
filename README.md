# scOntoMatch
scOntoMatch is an R package which unifies ontology annotation of scRNA-seq datasets to make them comparable across studies.


Maintainer: Yuyao Song ysong@ebi.ac.uk

Large-scale single cell atlases often have curated annotations using a standard ontology system. 
This package aims to match ontology labels between different datasets to make them comparable across studies. 

## Installation

```
## from source

install.packages("devtools")
devtools::install_github("YY-SONG0718/scOntoMatch")
library(scOntoMatch)
```

## Usage
First download the ontology `.obo` file from [the OBO foundry](https://obofoundry.org/). Common ontologies include:

  -  [Cell ontology](https://obofoundry.org/ontology/cl.html)
  -  [Uberon multi-species anatomy ontology](https://obofoundry.org/ontology/uberon.html) 
  -  [Drosophila gross anatomy (FBbt)](https://obofoundry.org/ontology/fbbt.html)
  -  [Zebrafish anatomy and development ontology (ZFA)](https://obofoundry.org/ontology/zfa.html)

The core functionality is to match ontology by direct mapping and mapping descendants to ancestor terms. 

```
ontoMatch(adata1, adata2, anno_col, onto_id_col, obo_file, propagate_relationships = c('is_a', 'part_of'))
```


This package also provides convenient plotting functions to help comprehend the hierarchy of cell types in any single cell dataset. 

```
plotOntoTree(ont, onts, plot_ancestors=TRUE, ont_query, fontsize = 20)
```

The [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) hosted at EBI provides uniformly analysed and annotated scRNA-Seq data across multiple species. 
Datasets with curated ontology labels are all great inputs to this package. scRNA-seq data stored as h5ad files can be downloaded via the [ftp site](http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/). 
Files that has extension `.project.h5ad` can be pass to `ontoMatch`. 


This package imports functions from [ontologyIndex](https://cran.r-project.org/web/packages/ontologyIndex/index.html), [ontologyPlot](https://cran.rstudio.com/web/packages/ontologyPlot/index.html) and 
[anndata](https://cran.r-project.org/web/packages/anndata/index.html)
