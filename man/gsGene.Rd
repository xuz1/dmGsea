\name{gsGene}
\alias{gsGene}
\title{
Gene set enrichment analysis based on gene level combined test P values
}
\description{
The function gsGene will first combine probe level P values to gene level P values, 
and then using threshold or ranking based methods to perform gene set enrichment analysis. 
If test data were provided together with test P values, correlations between probes will 
be adjusted in gene level P values. 
}
\usage{
gsGene(probe.p,Data4Cor=NULL,method="Threshold",FDRthre=0.05,nTopGene=NULL,
                   GeneProbeTable=NULL,arrayType=NULL,gSetName="KEGG",geneSet=NULL,
                   species="Human",combpMethod="fisher",combpAdjust = "nyholt",
                   outfile="gsGene",outGenep=FALSE,gseaParam=1,nperm=1e4,ncore=5)
}
\arguments{
  \item{probe.p}{A data frame that include columns for probe name "Name" and P values "p".} 
  \item{Data4Cor}{Methylation Data matrix (or an SummarizedExperiment object with 
	  assays(Data4Cor)$beta as the data matrix) used to generate probe.p, row as probe 
	  and column as samples. Row names should be the same with probe names in probe.p}
  \item{method}{Methods for testing gene set enrichment include options for "Threshold" and "Ranking".}
  \item{FDRthre}{False discovery rate threshold to select a list of significate genes for gene set 
  enrichment testing in threshold-based method.} 
  \item{nTopGene}{Specifies the number of top-ranked genes based on p-value. If provided, this will 
  override FDRthre argument to select gene list for gene set enrichment testing in threshold based method.}
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
  \item{outGenep}{TRUE or FALSE, if TRUE, combined gene level p values will be saved to "gsGene_genep.csv".}
  \item{gseaParam}{When method="Ranking", gene-level statistics -log(p) will be raised to the power of
          ‘gseaParam‘ to calculate enrichment scores}
  \item{nperm}{Number of permutation for Ranking method.}
  \item{ncore}{Number of cores will be used for computation}
}
\value{
	Results will be saved in files with name prefix specified by argument outfile.
}
\references{
Zongli Xu,Alison A. Motsinger-Reif, Liang Niu, Efficient gene set analysis for DNA methylation 
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
gsGene(probe.p=probe.p,Data4Cor=dat,GeneProbeTable=GeneProbeTable,
         method="Threshold",gSetName="KEGG",species="Human",ncore=50)
gsGene(probe.p=probe.p,Data4Cor=dat,GeneProbeTable=GeneProbeTable,
         method="Ranking",gSetName="KEGG",species="Human",outfile="geneRank",ncore=50)
gsGene(probe.p=probe.p,GeneProbeTable=GeneProbeTable,
         method="Threshold",gSetName="KEGG",species="Human",ncore=50)
}}


