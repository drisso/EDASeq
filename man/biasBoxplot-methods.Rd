\name{biasBoxplot-methods}
\docType{methods}
\alias{biasBoxplot}
\alias{biasBoxplot-methods}
\alias{biasBoxplot,numeric,numeric,numeric-method}
\alias{biasBoxplot,numeric,numeric-method}
\title{  Methods for Function \code{biasBoxplot} in Package \pkg{EDASeq} }
\description{
\code{biasBoxplot} produces a boxplot representing the distribution of a quantity of interest (e.g. gene counts, log-fold-changes, ...) stratified by a covariate (e.g. gene length, GC-contet, ...).
}

\usage{
biasBoxplot(x,y,num.bins,...)
}

\arguments{
\item{x}{A numeric vector with the quantity of interest (e.g. gene counts, log-fold-changes, ...)}
\item{y}{A numeric vector with the covariate of interest (e.g. gene length, GC-contet, ...)}
\item{num.bins}{A numeric value specifying the number of bins in wich to stratify \code{y}. Default to 10.}
\item{...}{See \code{\link{par}}}
}


\section{Methods}{
\describe{

\item{\code{signature(x = "numeric", y = "numeric", num.bins = "numeric")}}{
It plots a line representing the regression of every column of the matrix \code{x} on the numeric covariate \code{y}. One can pass the usual graphical parameters as additional arguments (see \code{\link{par}}).
}

}}
\keyword{methods}

\examples{
library(yeastRNASeq)
data(geneLevelData)
data(yeastGC)

sub <- intersect(rownames(geneLevelData), names(yeastGC))

mat <- as.matrix(geneLevelData[sub,])

data <- newSeqExpressionSet(mat,
                            phenoData=AnnotatedDataFrame(
                                      data.frame(conditions=factor(c("mut", "mut", "wt", "wt")),
                                                 row.names=colnames(geneLevelData))),
                            featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))

lfc <- log(geneLevelData[sub, 3] + 1) - log(geneLevelData[sub, 1] + 1)

biasBoxplot(lfc, yeastGC[sub], las=2, cex.axis=.7)

}

