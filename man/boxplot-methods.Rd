\name{boxplot-methods}
\docType{methods}
\alias{boxplot-methods}
\alias{boxplot,FastqQuality-method}
\alias{boxplot,SeqExpressionSet-method}
\title{  Methods for Function \code{boxplot} in Package \pkg{EDASeq} }
\description{
High-level functions to produce boxplots of some complex objects.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "FastqQuality")}}{
It plots the distribution of the quality per read position.
}

\item{\code{signature(x = "SeqExpressionSet")}}{
It plots the distribution of the log counts in each lane of \code{x}.
}
}}
\keyword{methods}

