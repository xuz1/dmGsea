\name{gsProbe}
\alias{gsProbe}
\title{
Gene set enrichment analysis based on probe level P values
}
\description{
The function gsProbe will not combine P values, instead, it will
adjust for the bias caused by variable number of probes per gene
using a noncentral hypergeometric test in gene set enrichment analysis.
}
\usage{
gsProbe(probe.p,FDRthre=0.05,nTopProbe=NULL,sigProbe=NULL,allProbe=NULL,
        GeneProbeTable=NULL,arrayType=NULL,gSetName=NULL,geneSet=NULL,
        species="Human",outfile="gsProbe",ncore=2)
}

\arguments{
    \item{probe.p}{A data frame that include columns for probe name "Name"
        and P values "p".} 
    \item{FDRthre}{False discovery rate threshold to select list of significate
        genes for gene set enrichment testing in threshold-based method.} 
    \item{nTopProbe}{Specifies the number of top-ranked probe based on p-value. 
        If provided, this will override FDRthre argument to select gene list
        for gene set enrichment testing.}
    \item{sigProbe}{A character vector to specifies a list significant probes,
        this will overrid FDRthre argument to select gene list for gene set 
        enrichment testing.}
    \item{allProbe}{Probe universe, use together with sigProbe argument.}
    \item{GeneProbeTable}{A data frame for probe to gene annotation with columns
        for probe names "Name" and entrez gene id "entrezid"}
    \item{arrayType}{Specifies array type to extract GeneProbeTable, currently 
        support DNA methylation array ("450K","EPIC"). User need to provide 
        GeneProbeTable for other types of dataset.}
    \item{gSetName}{Specifies gene set names, with options including "KEGG","GO"
        and "MSigDB", or subsets like "GO|BP".}
    \item{geneSet}{User-provided gene sets, where each set contains a list of 
        Entrez gene IDs, and the list name corresponds to the gene set name.}
    \item{species}{Specifies species including "Human" and "Mouse"}
    \item{outfile}{Prefix of output files.}
    \item{ncore}{Number of cores will be used for computation}
}
\value{
        Results will be saved to files with name prefix specified by outfile
        argument.
}
\references{
Zongli Xu,Alison A. Motsinger-Reif, Liang Niu, OmicGsea: Efficient Gene Set 
Enrichment Analysis for Complex Omics Datasets. in review 2024.
}

\author{Zongli Xu}

\examples{
if(FALSE){
kegg <- getKEGG(species="Human")
gene1 <- unique(as.vector(unlist(kegg[1:5])))
gene2 <- unique(as.vector(unlist(kegg[6:length(kegg)])))
gene1 <- rep(gene1,sample(1:10,length(gene1),replace=TRUE))
gene2 <- rep(gene2,sample(1:10,length(gene2),replace=TRUE))
p1 <- runif(length(gene1))*(1e-6)
p2 <- runif(length(gene2))
geneid <- c(gene1,gene2)
p <- c(p1,p2)
Name <- paste0("cg",1:length(p))
probe.p <- data.frame(Name=Name,p=p)
GeneProbeTable <- data.frame(Name=Name,entrezid=geneid)
dat <- matrix(runif(length(p)*100),ncol=100)
rownames(dat) <- Name
gsProbe(probe.p=probe.p,GeneProbeTable=GeneProbeTable,
    gSetName="KEGG",species="Human",ncore=2)
}}


