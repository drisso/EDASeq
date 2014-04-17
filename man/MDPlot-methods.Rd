\name{MDPlot-methods}
\docType{methods}
\alias{MDPlot}
\alias{MDPlot-methods}
\alias{MDPlot,matrix,numeric-method}
\alias{MDPlot,SeqExpressionSet,numeric-method}
\title{  Methods for Function \code{MDPlot} in Package \pkg{EDASeq} }
\description{
\code{MDPlot} produces a mean-difference smooth scatterplot of two lanes in an experiment.
}

\usage{
MDPlot(x,y,...)
}

\arguments{
\item{x}{Either a numeric matrix or a \code{\linkS4class{SeqExpressionSet}} object containing the gene expression.
}
\item{y}{A numeric vecor specifying the lanes to be compared.}
\item{...}{See \code{\link{par}}}
}

\details{
The mean-difference (MD) plot is a useful plot to visualize difference in two lanes of an experiment. From a MDPlot one can see if normalization is needed and if a linear scaling is sufficient or nonlinear normalization is more effective.

The MDPlot also plots a lowess fit (in red) underlying a possible trend in the bias related to the mean expression.

}
\section{Methods}{
\describe{

\item{\code{signature(x = "matrix", y = "numeric")}}{
}

\item{\code{signature(x = "SeqExpressionSet", y = "numeric")}}{

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

MDPlot(data,c(1,3))

}

