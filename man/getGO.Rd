\name{getGO}
\alias{getGO}
\title{
Extract GO pathway gene set
}
\description{
    Extract GO pathway gene set.
}
\usage{
getGO(subset="BP",species="Human")
}

\arguments{
  \item{species}{Specifies species including "Human" and "Mouse"}
  \item{subset}{Specifies subset of GO pathway including "BP","MF" and "CC"}
}
\value{
    A list of pathways and entrezid id in each pathway.
}

\author{
Zongli Xu
}
\examples{
GO <- getGO(subset="BP",species="Human")
}


