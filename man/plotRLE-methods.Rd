\name{plotRLE-methods}
\docType{methods}
\alias{plotRLE}
\alias{plotRLE-methods}
\alias{plotRLE,matrix-method}
\alias{plotRLE,SeqExpressionSet-method}
\title{  Methods for Function \code{plotRLE} in Package \pkg{EDASeq} }
\description{
\code{plotRLE} produces a Relative Log Expression (RLE) plot of the counts in \code{x}
}

\usage{
plotRLE(x, ...)
}

\arguments{
\item{x}{Either a numeric matrix or a \code{\linkS4class{SeqExpressionSet}} object containing the gene expression.
}
\item{...}{See \code{\link{par}}}
}

\details{
The Relative Log Expression (RLE) plot is a useful diagnostic plot to visualize the differences between the distributions
of read counts across samples. 

It shows the boxplots of the log-ratios of the gene-level read counts of each sample to those of a reference sample (defined as the median across the samples). Ideally, the distributions should be centered around the zero line and as tight as possible. Clear deviations indicate the need for normalization and/or the presence of outlying samples.
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


plotRLE(data, col=rep(2:3, each=2))

}

