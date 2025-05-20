\name{getReactome}
\alias{getReactome}
\title{
Extract Reactome pathway gene set
}
\description{
    Extract Reactome pathway gene set.
}
\usage{
getReactome(species="Human")
}

\arguments{
    \item{species}{Specifies species including "Human" and "Mouse"}
}
\value{
    A list of pathways and entrezid id in each pathway.
}
\references{
Zongli Xu, Liang Niu, OmicGsea: Efficient gene set analysis for 
complex Omics data. sutmitted 2024.
}

\author{Zongli Xu}

\examples{
Reactome <- getReactome(species="Human")
}


