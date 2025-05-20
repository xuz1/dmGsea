\name{getMSigDB}
\alias{getMSigDB}
\title{
Extract MsigDB pathway gene set
}
\description{
    Extract MsigDB pathway gene set.
}
\usage{
getMSigDB(subset="C2",species="Human")
}

\arguments{
    \item{species}{Specifies species including "Human" and "Mouse"}
    \item{subset}{Specifies subset of MsigDB pathway including "H",
    "C1","C2","C3","C4","C5","C6","C7" and"C8"}
}
\value{
    A list of pathways and entrezid id in each pathway.
}

\author{Zongli Xu}

\examples{
MsigDB <- getMSigDB(subset="C2",species="Human")
}


