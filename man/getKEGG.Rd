\name{getKEGG}
\alias{getKEGG}
\title{
Extract KEGG pathway gene set
}
\description{
	Extract KEGG pathway gene set.
}
\usage{
getKEGG(species="Human")
}

\arguments{
  \item{species}{Specifies species including "Human" and "Mouse"}
}
\value{
	A list of pathways and entrezid id in each pathway.
}

\author{Zongli Xu}

\examples{
kegg=getKEGG(species="Human")
}


