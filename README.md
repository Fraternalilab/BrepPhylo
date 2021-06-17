# BrepPhylo

An R package with functions to perform phylogenetic reconstruction of B-cell receptor (BCR) lineages.

The `BrepPhylo` package contains a set of functions written for the purpose of performing lineage reconstruction of a large number of B cell clonotypes obtained via high-throughput sequencing of immunoglobulin repertoires. The aim is to provide functionalities as 'wrappers' around popular, existing tools for lineage reconstruction, provide functionalities to detect Class-Switch Recombination (CSR) events from lineage trees, and provide easy-to-interpret output to analyse these lineages.  

## Installation

In R:

```
require(devtools)
install_github("Fraternalilab/BrepPhylo", ref = "main", dependencies = TRUE)
```

## Examples

Vignettes in [HTML](http://htmlpreview.github.io/?https://github.com/Fraternalilab/BrepPhylo/blob/main/vignettes/phylogenetic.html) and [PDF](https://github.com/Fraternalilab/BrepPhylo/blob/main/vignettes/phylogenetic.pdf) formats are provided here.
