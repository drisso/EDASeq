\name{biasPlot-methods}
\docType{methods}
\alias{biasPlot}
\alias{biasPlot-methods}
\alias{biasPlot,matrix,numeric-method}
\alias{biasPlot,SeqExpressionSet,character-method}
\title{  Methods for Function \code{biasPlot} in Package \pkg{EDASeq} }
\description{
\code{biasPlot} produces a plot of the \code{\link{lowess}} regression of the counts on a covariate of interest, tipically the GC-content or the length of the genes.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "matrix", y = "numeric")}}{
It plots a line representing the regression of every column of the matrix \code{x} on the numeric covariate \code{y}. One can pass the usual graphical parameters as additional arguments (see \code{\link{par}}).
}

\item{\code{signature(x = "SeqExpressionSet", y = "character")}}{
It plots a line representing the regression of every lane in \code{x} on the covariate specified by \code{y}. \code{y} must be one of the column of the \code{featureData} slot of the \code{x} object. One can pass the usual graphical parameters as additional arguments (see \code{\link{par}}). The parameter \code{color_code} (optional) must be a number specifying the column of \code{phenoData} to be used for color-coding. By default it is color-coded according to the first column of \code{phenoData}. If \code{legend=TRUE} and \code{col} is not specified a legend with the information stored in \code{phenoData} is added.
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

biasPlot(data,"gc",ylim=c(0,5),log=TRUE)

}

