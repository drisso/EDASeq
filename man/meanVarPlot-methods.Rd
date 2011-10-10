\name{meanVarPlot-methods}
\docType{methods}
\alias{meanVarPlot}
\alias{meanVarPlot-methods}
\alias{meanVarPlot,SeqExpressionSet-method}
\title{ Methods for Function \code{meanVarPlot} in Package \pkg{EDASeq}}
\description{
\code{meanVarPlot} produces a \code{smoothScatter} plot of the mean variance relation.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "SeqExpressionSet")}}{
It takes as additional argument \code{log}, which if true consider the logarithm of the counts before computing mean and variance. To avoid missing values, we consider the maximum between 0 and the log of the counts. Along with the scatter plot the function plots a line representing the \code{\link{lowess}} fit.
}
}}
\keyword{methods}

