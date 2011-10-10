\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,BamFileList,FastqFileList-method}
\title{  Methods for Function \code{plot} in Package \pkg{EDASeq} }
\description{
High-level function to produce plots given one \code{BamFileList} object and one \code{FastqFileList} object.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "BamFileList", y = "FastqFileList")}}{
It produce a barplot of the percentage of mapped reads. If \code{strata=TRUE} it stratifies the bars according to the unique/non-unique mapped reads.
To be meaningful, \code{x} should be a set of aligned reads and \code{y} a set of raw reads on the same samples.
}
}}
\keyword{methods}
