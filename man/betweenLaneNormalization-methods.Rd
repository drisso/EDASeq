\name{betweenLaneNormalization-methods}
\docType{methods}
\alias{betweenLaneNormalization}
\alias{betweenLaneNormalization-methods}
\alias{betweenLaneNormalization,matrix-method}
\alias{betweenLaneNormalization,SeqExpressionSet-method}
\title{ Methods for Function \code{betweenLaneNormalization} in Package \pkg{EDASeq} }
\description{
 Between-lane normalization for sequencing depth and possibly other distributional differences between lanes.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "matrix")}}{
It returns a matrix with the normalized counts if \code{offset=FALSE} or with the offset if \code{offset=TRUE}.
}

\item{\code{signature(x = "SeqExpressionSet")}}{
It returns a \code{linkS4class{SeqExpressionSet}} with the normalized counts in the \code{normalizedCounts} slot and with the offset in the \code{offset} slot (if \code{offset=TRUE}).
}
}}

\usage{
betweenLaneNormalization(x, which=c("median","upper","full"), offset=FALSE, round=TRUE)
}

\arguments{
\item{x}{A numeric matrix representing the counts or a \code{\linkS4class{SeqExpressionSet}} object.}
\item{which}{Method used to normalized. See the details section and the reference below for details.}
\item{offset}{Should the normalized value be returned as an offset leaving the original counts unchanged?}
\item{round}{If TRUE the normalization returns rounded values (pseudo-counts). Ignored if offset=TRUE.}
}



\details{
This method implements three normalizations described in Bullard et al. (2010).
The methods are:
\describe{
\item{\code{median}:}{a scaling normalization that forces the median of each lane to be the same.}
\item{\code{upper}:}{the same but with the upper quartile.}
\item{\code{full}:}{a non linear full quantile normalization, in the spirit of the one used in microarrays.}
}
}

\author{
Davide Risso.
}

\references{
J. H. Bullard, E. A. Purdom, K. D. Hansen and S. Dudoit (2010). Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics Vol. 11, Article 94.

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

norm <- betweenLaneNormalization(data, which="full", offset=FALSE)

}

\keyword{methods}

