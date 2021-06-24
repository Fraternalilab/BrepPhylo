# BrepPhylo

An R package with functions to perform phylogenetic reconstruction of B-cell receptor (BCR) lineages.

The `BrepPhylo` package contains a set of functions written for the purpose of performing lineage reconstruction of a large number of B cell clonotypes obtained via high-throughput sequencing of immunoglobulin repertoires. The aim is to provide 'wrappers' around popular, existing tools for lineage reconstruction, detect Class-Switch Recombination (CSR) events from lineage trees, and generate easy-to-interpret output to analyse these lineages.  

## Installation

In R:

```
require(devtools)
install_github("Fraternalilab/BrepPhylo", ref = "main", dependencies = TRUE)
```

To run the examples in the vignettes the `dnapars` program from the [`PHYLIP`](https://evolution.genetics.washington.edu/phylip.html) package is required. An executable of the program is shipped with the package; **this works for Linux systems**. For other operating systems please follow instructions from PHYLIP to download and install the program.

## Examples

Vignettes in [HTML](http://htmlpreview.github.io/?https://github.com/Fraternalilab/BrepPhylo/blob/main/vignettes/phylogenetic.html) and [PDF](https://github.com/Fraternalilab/BrepPhylo/blob/main/vignettes/phylogenetic.pdf) formats are provided here.

## Contact

Joseph Ng ([@josef0731](https://github.com/josef0731)). For issues with this package please open an issue in this github repository.

## Acknowledgement

To members of the Dunn-Walters lab, University of Surrey for developing ideas to implement in the package, and to Christian Margreitter ([@CMargreitter](https://github.com/CMargreitter)) for building the initial version of the package.
