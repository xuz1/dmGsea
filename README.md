# dmGsea
## Flexible Gene Set Enrichment Analysis for DNA Methylation Data
## Introduction
The R package <b>dmGsea</b> provides gene set enrichment analysis functions for complex Omics data. In many omics datasets, such as Illumina DNA methylation data, genes are often represented by variable numbers of correlated probes, and a single probe may map to multiple genes. This complex data structure poses significant challenges for gene set enrichment analysis. The dmGsea package addresses these challenges by providing both threshold- and ranking-based competitive enrichment analysis methods, effectively accounting for probe correlation and probe number bias. These versatile functions can be applied to a wide range of omics datasets, including DNA methylation, proteomics and gene expression data.

## Functions

<ul>
<li>`gsGene()`:  Gene set enrichment analysis based on gene level combined p-values</li>
<li>`gsPG()`:  Gene set enrichment analysis based on combined probe group p-values</li>
<li>`gsProbe()`:  Gene set enrichment analysis based on probe level p-values</li>
<li>`gsRank()`:  Ranking-based gene set enrichment analysis based on gene level p-values</li>
</ul>
