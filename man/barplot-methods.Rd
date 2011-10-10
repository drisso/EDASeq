\name{barplot-methods}
\docType{methods}
\alias{barplot-methods}
\alias{barplot,BamFile-method}
\alias{barplot,BamFileList-method}
\alias{barplot,FastqFileList-method}
\title{  Methods for Function \code{barplot} in Package \pkg{EDASeq} }
\description{
High-level functions to produce barplots of some complex objects.
}
\section{Methods}{
\describe{

\item{\code{signature(height = "BamFile")}}{
Usage: barplot(height,strata=c("rname","strand"))
It produces a barplot of the total number of reads in each chromosome (if "rname") or strand.
}

\item{\code{signature(height = "BamFileList")}}{
It produces a barplot of the total number of reads in each object in \code{height}. If \code{unique=TRUE} is specified, it stratified the total by uniquely/non-uniquely mapped reads.
}

\item{\code{signature(height = "FastqFileList")}}{
It produces a barplot of the total number of reads in each object in \code{height}.
}

}}
\keyword{methods}

