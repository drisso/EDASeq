\name{plotPCA-methods}
\docType{methods}
\alias{plotPCA}
\alias{plotPCA-methods}
\alias{plotPCA,matrix-method}
\alias{plotPCA,SeqExpressionSet-method}
\title{  Methods for Function \code{plotPCA} in Package \pkg{EDASeq} }
\description{
\code{plotPCA} produces a Principal Component Analysis (PCA) plot of the counts in \code{object}
}

\usage{
\S4method{plotPCA}{matrix}(object, k=2, labels=TRUE, isLog=FALSE, ...)
\S4method{plotPCA}{SeqExpressionSet}(object, k=2, labels=TRUE, ...)
}

\arguments{
\item{object}{Either a numeric matrix or a \code{\linkS4class{SeqExpressionSet}} object containing the gene expression.}

\item{k}{The number of principal components to be plotted.}

\item{labels}{Logical. If \code{TRUE}, and \code{k=2}, it plots the \code{colnames} of \code{object} as point labels.}

\item{isLog}{Logical. Set to \code{TRUE} if the data are already on the log scale.}

\item{...}{See \code{\link{par}}}
}

\details{
The Principal Component Analysis (PCA) plot is a useful diagnostic plot to highlight differences in the distribution of replicate samples, by projecting the samples into a lower dimensional space.

If there is strong differential expression between two classes, one expects the samples to cluster by class in the first few Principal Components (PCs) (usually 2 or 3 components are enough). This plot also highlights possible batch effects and/or outlying samples.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "matrix")}}{
}

\item{\code{signature(x = "SeqExpressionSet")}}{

}
}}
\keyword{methods}

\examples{
library(yeastRNASeq)
data(geneLevelData)

mat <- as.matrix(geneLevelData)

data <- newSeqExpressionSet(mat,
                            phenoData=AnnotatedDataFrame(
                                      data.frame(conditions=factor(c("mut", "mut", "wt", "wt")),
                                                 row.names=colnames(geneLevelData))))

plotPCA(data, col=rep(1:2, each=2))

}

