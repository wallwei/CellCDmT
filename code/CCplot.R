# install development version
library('devtools')

devtools::install_github("Sarah145/CCPlotR")

# or install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("CCPlotR")
