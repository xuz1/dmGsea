\name{gsPG}
\alias{gsPG}
\title{
Gene set enrichment analysis based on combined probe group P values
}
\description{
The function gsPG will first group probes with same gene annotation and combine P values 
for each group, and then using noncentral hypergeometric or Monte Carlo method to perform 
gene set enrichment analysis. If test data were provided together with test P values, 
correlations between probes will be adjusted in probe group level P values. 
}
\usage{
gsPG(probe.p,Data4Cor=NULL,FDRthre=0.05,nTopPG=NULL,MonteCarlo=FALSE,
                 GeneProbeTable=NULL,arrayType=NULL,gSetName=NULL,geneSet=NULL,
                 species="Human",combpMethod="fisher",combpAdjust = "nyholt",
                 outfile="gsPG",ncore=5)
}

\arguments{
  \item{probe.p}{A data frame that include columns for probe name "Name" and P values "p".} 
  \item{Data4Cor}{Methylation Data matrix (or an SummarizedExperiment object with 
	  assays(Data4Cor)$beta as the data matrix)used to generate probe.p, row as probe
	  and column as samples. Row names should be the same with probe names in probe.p}
  \item{FDRthre}{False discovery rate threshold to select a list of significate genes for gene set 
  enrichment testing in threshold-based method.} 
  \item{nTopPG}{Specifies the number of top-ranked probe group based on p-value. If provided, this will 
  override FDRthre argument to select gene list for gene set enrichment testing in threshold based method.}
  \item{MonteCarlo}{If TRUE, a Monte Carlo method will be used to perform geneset test, otherwise, a 
	  noncentral hypergeometric test will be used.} 
  \item{GeneProbeTable}{A data frame for probe to gene annotation with columns for probe names "Name" and 
  entrez gene id "entrezid"}
  \item{arrayType}{Specifies array type to extract GeneProbeTable, currently support DNA methylation array
  ("450K","EPIC"). User need to provide GeneProbeTable for other types of dataset.}
  \item{gSetName}{Specifies gene set names, with options including "KEGG","GO" and "MSigDB", 
  or subsets like "GO|BP".}
  \item{geneSet}{User-provided gene sets, where each set contains a list of Entrez gene IDs, and the 
  list name corresponds to the gene set name.}
  \item{species}{Specifies species including "Human" and "Mouse"}
  \item{combpMethod}{Specifies combine p-value methods including "fisher", "invchisq", "stouffer" and "tippett"}
  \item{combpAdjust}{Specifies method to adjust for dependence between probes including "none", 
  "nyholt", "liji","gao" and "galwey".}
  \item{outfile}{Prefix of output files.}
  \item{ncore}{Number of cores will be used for computation}
}
\value{
	Results will be saved to files with name prefix specified by outfile argument.
}
\references{
Zongli Xu, Alison A. Motsinger-Reif, Liang Niu, Efficient gene set analysis for DNA methylation 
addressing probe dependency and bias. in review 2024.
}

\author{Zongli Xu}

\examples{
if(FALSE){
kegg=getKEGG(species="Human")
gene1=unique(as.vector(unlist(kegg[1:5])))
gene2=unique(as.vector(unlist(kegg[6:length(kegg)])))
gene1=rep(gene1,sample(1:10,length(gene1),replace=TRUE))
gene2=rep(gene2,sample(1:10,length(gene2),replace=TRUE))
p1=runif(length(gene1))*(1e-3)
p2=runif(length(gene2))
geneid=c(gene1,gene2)
p=c(p1,p2)
Name=paste0("cg",1:length(p))
probe.p=data.frame(Name=Name,p=p)
GeneProbeTable=data.frame(Name=Name,entrezid=geneid)
dat=matrix(runif(length(p)*100),ncol=100)
rownames(dat)=Name
gsPG(probe.p=probe.p,Data4Cor=dat,GeneProbeTable=GeneProbeTable,
        gSetName="KEGG",species="Human",ncore=50)
gsPG(probe.p=probe.p,Data4Cor=dat,GeneProbeTable=GeneProbeTable,
         MonteCarlo=TRUE,gSetName="KEGG",species="Human",outfile="genePG_MC",ncore=50)
}}


