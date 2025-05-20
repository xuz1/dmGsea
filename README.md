# dmGsea
## Efficient Gene Set Enrichment Analysis for DNA Methylation Data
## Introduction
The R package <b>dmGsea</b> is specifically designed for DNA methylation data,
 provides functions to perform gene set enrichment analysis while addressing 
probe dependency and probe number bias. The package supports Illumina 450K, 
EPIC, and mouse methylation arrays and can be extended to other omics data 
with user-provided probe-to-gene mapping annotations.

## Installation

### From Bioconductor (recommended)
You can install the stable version of **dmGsea** from 
[Bioconductor](https://bioconductor.org) using the following command:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dmGsea")

