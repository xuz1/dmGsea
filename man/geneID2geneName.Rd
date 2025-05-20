\name{geneID2geneName}
\alias{geneID2geneName}
\title{
Extract gene names for given entrezids
}
\description{
    Extract gene name for given entrezid.
}
\usage{
geneID2geneName(gid="NULL",species="Human")
}

\arguments{
  \item{gid}{Specifies a vector of entrezids}
  \item{species}{Specifies species including "Human" and "Mouse"}
}
\value{
    A list of gene names.
}
\author{
Zongli Xu
}
\examples{
gname <- geneID2geneName(gid=c("672","675"),species="Human")
}


