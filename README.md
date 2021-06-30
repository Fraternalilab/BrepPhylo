# BrepPhylo

An R package with functions to perform phylogenetic reconstruction of B-cell receptor (BCR) lineages.

The `BrepPhylo` package contains a set of functions written for the purpose of performing lineage reconstruction of a large number of B cell clonotypes obtained via high-throughput sequencing of immunoglobulin repertoires. The aim is to provide 'wrappers' around popular, existing tools for lineage reconstruction, detect Class-Switch Recombination (CSR) events from lineage trees, and generate easy-to-interpret output to analyse these lineages.  

## Installation

### MacOS / Linux

In R:

```
# install these dependencies on Bioconductor if you haven't had them installed already
install.packages('BiocManager')
BiocManager::install( c("Biostrings", "msa", "GenomicRanges", "GenomicAlignments") )

# install the package; the following install those dependent packages available on CRAN
require(devtools)
install_github("Fraternalilab/BrepPhylo", ref = "main", dependencies = TRUE)
```

To run the examples in the vignettes the `dnapars` program from the [`PHYLIP`](https://evolution.genetics.washington.edu/phylip.html) package is required. An executable of the program is shipped with the package; **this works for Linux systems**. 

For MacOS please follow instructions from PHYLIP to download and install the program. You will need to supply the **absolute** filepath of the `dnapars` executable when you run some of the functions in `BrepPhylo`. See the vignettes (below) for details.

### Windows

1. Make sure you have Rtools installed. Instructions can be found [here](https://cran.r-project.org/bin/windows/Rtools/). Basically:
  + Run the installer. You can install the program anywhere you want.
  + from R, run: `writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")`
  + **Close and restart R**. Verify this is set up by running this in R `Sys.which("make")`. It should return the path where Rtools is installed.
2. Install the `devtools` R package. 
  + In R: `install.packages("devtools")`
3. Next download PHYLIP from the link above. Make sure you selected the Windows pre-computed binaries. 
  + You will then download a zip file. Extract the zip file.
  + The file `dnapars.exe` in the `exe` subfolder of the extracted directory is the filepath you need to supply when running certain functions in `BrepPhylo`. See the vignettes (below) for details.
4. `BrepPhylo` calls the `alakazam` R package in running PHYLIP. As of 30th June 2021, you need to download the development version of the alakazam package rather than the standard avaialble pre-compiled version (**only for Windows; the CRAN available version works as expected for MacOS/Linux machines**).
  + In R: `devtools::install_bitbucket("kleinstein/alakazam")`
5. (*Optional*) If you wish to generate alignment PDF files for sequences, you need a working copy of LaTex. Check out distributions such as [MiKTeX](http://miktex.org/) or [TeX Live](http://www.tug.org/texlive) and install them. **Note:** this is an optional step, you can indicate not to generate such output when running `BrepPhylo` functions, lineage trees will still be reconstructed and you can still use other functionalities in the package as designed.
6. Now follow the instructions above for MacOS/Linux to install `BrepPhylo`.

We have tested running the vignette on Windows machines and it works as expected.

## Examples

Vignettes in [HTML](http://htmlpreview.github.io/?https://github.com/Fraternalilab/BrepPhylo/blob/main/vignettes/phylogenetic.html) and [PDF](https://github.com/Fraternalilab/BrepPhylo/blob/main/vignettes/phylogenetic.pdf) formats are provided here.

## Contact

Joseph Ng ([@josef0731](https://github.com/josef0731)). For issues with this package please open an issue in this github repository.

## Acknowledgement

To members of the Dunn-Walters lab, University of Surrey for developing ideas to implement in the package, and to Christian Margreitter ([@CMargreitter](https://github.com/CMargreitter)) for building the initial version of the package.
