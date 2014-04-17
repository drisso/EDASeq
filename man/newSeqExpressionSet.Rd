\name{newSeqExpressionSet}
\alias{newSeqExpressionSet}

\title{
Function to create a new \code{\linkS4class{SeqExpressionSet}} object.
}
\description{
User-level function to create new objects of the class \code{\linkS4class{SeqExpressionSet}}.
}
\usage{
newSeqExpressionSet(counts,
                    normalizedCounts = matrix(data=NA, nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts)),
                    offset = matrix(data=0, nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts)),
                    phenoData = annotatedDataFrameFrom(counts, FALSE),
                    featureData = annotatedDataFrameFrom(counts, TRUE),
                    ...)
}
\arguments{
  \item{counts}{
  A matrix containing the counts for an RNA-Seq experiment. One column for each lane and one row for each gene.
}
  \item{normalizedCounts}{
   A matrix with the same dimensions of \code{counts} with the normalized counts.
}
  \item{offset}{
 A matrix with the same dimensions of \code{counts} defining the offset (usually useful for normalization purposes). See the package vignette for a discussion on the offset.
}
  \item{phenoData}{
A data.frame or \code{\link{AnnotatedDataFrame}} with sample information, such as biological condition, library preparation protocol, flow-cell,...
}
  \item{featureData}{
A data.frame or \code{\link{AnnotatedDataFrame}} with feature information, such as gene length, GC-content, ...
}
  \item{\dots}{
 Other arguments will be passed to the constructor inherited from \code{\linkS4class{eSet}}.
}
}

\value{
An object of class \code{\linkS4class{SeqExpressionSet}}.
}

\author{
Davide Risso
}

\seealso{
\code{\linkS4class{SeqExpressionSet}}
}
\examples{

counts <- matrix(data=0, nrow=100, ncol=4)
for(i in 1:4) {
counts[, i] <- rpois(100, lambda=50)
}
cond <- c(rep("A", 2), rep("B", 2))

counts <- newSeqExpressionSet(counts, phenoData=data.frame(conditions=cond))

}
