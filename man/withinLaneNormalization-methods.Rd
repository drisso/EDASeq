\name{withinLaneNormalization-methods}
\docType{methods}
\alias{withinLaneNormalization}
\alias{withinLaneNormalization-methods}
\alias{withinLaneNormalization,matrix,numeric-method}
\alias{withinLaneNormalization,SeqExpressionSet,character-method}
\title{  Methods for Function \code{withinLaneNormalization} in Package \pkg{EDASeq} }
\description{
 Within-lane normalization for GC-content (or other lane-specific) bias.
}

\section{Methods}{
\describe{
\item{\code{signature(x = "matrix", y = "numeric")}}{

It returns a matrix with the normalized counts if \code{offset=FALSE} or with the offset if \code{offset=TRUE}.

}

\item{\code{signature(x = "SeqExpressionSet", y = "character")}}{

It returns a \code{\linkS4class{SeqExpressionSet}} with the normalized counts in the \code{normalizedCounts} slot and with the offset in the \code{offset} slot (if \code{offset=TRUE}).

}
}}

\usage{
withinLaneNormalization(x, y, which=c("loess","median","upper","full"), offset=FALSE, num.bins=10, round=TRUE)
}

\arguments{
\item{x}{A numeric matrix representing the counts or a \code{\linkS4class{SeqExpressionSet}} object.}
\item{y}{A numeric vector representing the covariate to normalize for (if \code{x} is a matrix) or a character vector with the name of the covariate (if \code{x} is a \code{\linkS4class{SeqExpressionSet}} object). Usually it is the GC-content.}
\item{which}{Method used to normalized. See the details section and the reference below for details.}
\item{offset}{Should the normalized value be returned as an offset leaving the original counts unchanged?}
\item{num.bins}{The number of bins used to stratify the covariate for \code{median}, \code{upper} and \code{full} methods. Ignored if \code{loess}. See the reference for a discussion on the number of bins.}
\item{round}{If TRUE the normalization returns rounded values (pseudo-counts). Ignored if offset=TRUE.}
}

\details{
This method implements four normalizations described in Risso et al. (2011).

The \code{loess} normalization transforms the data by regressing the counts on \code{y} and subtracting the loess fit from the counts to remove the dependence.

The \code{median}, \code{upper} and \code{full} normalizations are based on the stratification of the genes based on \code{y}. Once the genes are stratified in \code{num.bins} strata, the methods work as follows. 
\describe{
\item{\code{median}:}{scales the data to have the same median in each bin.}
\item{\code{upper}:}{the same but with the upper quartile.}
\item{\code{full}:}{forces the distribution of each stratum to be the same using a non linear full quantile normalization, in the spirit of the one used in microarrays.}
}
}

\author{
Davide Risso.
}

\references{
D. Risso, K. Schwartz, G. Sherlock and S. Dudoit (2011). GC-Content Normalization for RNA-Seq Data. Manuscript in Preparation.
}

\examples{
library(yeastRNASeq)
data(geneLevelData)
data(yeastGC)

sub <- intersect(rownames(geneLevelData), names(yeastGC))

mat <- as.matrix(geneLevelData[sub, ])

data <- newSeqExpressionSet(mat,
                            phenoData=AnnotatedDataFrame(
                                      data.frame(conditions=factor(c("mut", "mut", "wt", "wt")),
                                                 row.names=colnames(geneLevelData))),
                            featureData=AnnotatedDataFrame(data.frame(gc=yeastGC[sub])))

norm <- withinLaneNormalization(data, "gc", which="full", offset=FALSE)

}
\keyword{methods}

