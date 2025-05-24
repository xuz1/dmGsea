\name{gsRank}
\alias{gsRank}
\title{
Gene set enrichment analysis based on gene level P values
}
\description{
Similar to GSEA and fgsea, the function gsRank will perform geneset enrichment
analysis based on pre-ranked gene list (such as results from gene expression
array) using ranking based method.
}
\usage{
gsRank(stats,outfile="gsRank",scoreType="std",gSetName=NULL,
        geneSet=NULL,gseaParam=1,species="Human",nperm=1e4,ncore=1)
}

\arguments{
    \item{stats}{Named vector of gene-level test statistics, Names should be the
        same with names in geneSet.} 
    \item{outfile}{Prefix of output files.}
    \item{scoreType}{Defines the enrichment test type, two-sided as in the 
        original GSEA ("std"), or one-tailed tests "pos" or "neg".}
    \item{gSetName}{Specifies gene set names, with options including "KEGG",
        "GO" and "MSigDB", or subsets like "GO|BP".}
    \item{geneSet}{User-provided gene sets, where each set contains a list of 
        Entrez gene IDs, and the list name corresponds to the gene set name.}
    \item{gseaParam}{gene-level statistics will be raised to the power of 
        ‘gseaParam‘ to calculate enrichment scores}
    \item{species}{Specifies species including "Human" and "Mouse"}
    \item{nperm}{Specifies the number of permutations}
    \item{ncore}{Number of cores will be used for computation}
}
\value{
        Results will be saved to files with name prefix specified by outfile
            argument.
}
\references{
Zongli Xu, Alison A. Motsinger-Reif, Liang Niu, OmicGsea: Efficient Gene Set 
Enrichment Analysis for Complex Omics Datasets. in review 2024.
}
\author{Zongli Xu}

\examples{
if(FALSE){
kegg <- getKEGG(species="Human")
gene <- unique(as.vector(unlist(kegg)))
p <- runif(length(gene))
names(p) <- gene
stats <- -log(p)*sample(c(1,-1),length(p),replace=TRUE)

#traditional GSEA analysis, enrichment toward higher or lower end of statstics
stats <- sort(stats,decr=TRUE)
gsRank(stats=stats,gSetName="KEGG",scoreType="std",outfile="gsea9",nperm=1e5,
    ncore=1)
#enrichment of genes with higher statistics
stats <- sort(abs(stats),decr=TRUE)
gsRank(stats=stats,gSetName="KEGG",scoreType="std",outfile="gsea10",nperm=1e5,
    ncore=1)
}}
