\name{getIlluminaAnnotation}
\alias{getIlluminaAnnotation}
\title{
Extract CpG-Gene corresponding table from Illuminal annotation file.
}
\description{
Extract CpG-Gene correspondence table from Illuminal annotation file.
}
\usage{
getIlluminaAnnotation(arrayType=c("450K","EPIC"))
}

\arguments{
  \item{arrayType}{Specifies Illuminal methylation array type, including "450K" and "EPIC"}
}
\value{
	A data from table for CpG-Gene correspondence relationship.
}

\author{
Zongli Xu
}
\examples{
anno=getIlluminaAnnotation(arrayType="EPIC")
}


