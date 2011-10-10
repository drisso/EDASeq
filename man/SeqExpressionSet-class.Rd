\name{SeqExpressionSet-class}
\docType{class}
\alias{SeqExpressionSet-class}
\alias{offst}
\alias{offst<-}
\alias{initialize,SeqExpressionSet-method}
\alias{exprs,SeqExpressionSet-method}
\alias{exprs<-,SeqExpressionSet,ANY-method}
\alias{offst,SeqExpressionSet-method}
\alias{offst<-,SeqExpressionSet,ANY-method}
\alias{offst<-,SeqExpressionSet-method}
\alias{coerce,SeqExpressionSet,CountDataSet-method}



\title{"SeqExpressionSet" class for collections of short reads}

\description{

This class represents a collection of digital expression data (usually counts from RNA-Seq technology) along with sample information.

}

\section{Objects from the Class}{

  Objects of this class can be created from a call to the
  \code{\link{newSeqExpressionSet}} constructor.

}

\section{Extends}{
	Class \code{eSet}, directly.
	Class \code{VersionedBiobase}, by class \code{eSet}, distance 2.
	Class \code{Versioned}, by class \code{eSet}, distance 3.
}

\section{Slots}{
  Inherited from \code{eSet}:
  
  \describe{
     \item{\code{assayData}}{Contains matrices with equal dimensions, and with
          column number equal to \code{nrow(phenoData)}.\code{assayData} must
          contain a matrix \code{exprs} with rows represening features
          (e.g., genes) and columns representing samples.
          The optional matrix \code{offset} can be added to represent a normalization offset to be used for differential expression analysis. See the vignette for details.
          Class: \code{\link{AssayData-class}}.}

    \item{\code{phenoData}}{Sample information. For compatibility with DESeq, there should be at least the column \code{conditions}. See \code{\link{eSet}} for details.}

    \item{\code{featureData}}{Feature information. It is recomended to include at least length and GC-content information. This slot is used for \code{\link{withinLaneNormalization}}. See \code{\link{eSet}} for details.}

    \item{\code{experimentData}}{See \code{\link{eSet}}}

    \item{\code{annotation}}{See \code{\link{eSet}}}

    \item{\code{protocolData}}{See \code{link{eSet}}}
    }
}


\section{Methods}{
  See \code{\link{eSet}} for inherited methods. Additional methods:

  \describe{

    \item{exprs}{\code{signature(object="SeqExpressionSet")}: returns the \code{exprs} matrix. }

    \item{exprs<-}{\code{signature(object = "SeqExpressionSet")}: method to replace the \code{exprs} matrix.}

    \item{offst}{\code{signature(object = "SeqExpressionSet")}:  returns the \code{offset} matrix.}

    \item{offst<-}{\code{signature(object = "SeqExpressionSet")}: method to replace the \code{offset} slot.}

    \item{boxplot}{\code{signature(x = "SeqExpressionSet")}: produces a boxplot of the log counts.}

    \item{meanVarPlot}{\code{signature(x = "SeqExpressionSet")}: produces a \code{\link{smoothScatter}} plot of the mean variance relation. See \code{\link{meanVarPlot}} for details.}

    \item{biasPlot}{\code{signature(x = "SeqExpressionSet", y = "character")}: produces a plot of the \code{\link{lowess}} regression of the counts on some covariate of interest (usually GC-content or length). See \code{\link{biasPlot}} for details.}

    \item{wihtinLaneNormalization}{\code{signature(x = "SeqExpressionSet", y = "missing")}: within lane normalization for GC-content (or other lane specific) bias. See \code{\link{withinLaneNormalization}} for details.}

    \item{betweenLaneNormalization}{\code{signature(x = "SeqExpressionSet")}: between lane normalization for sequencing depth and possibly other distributional differences between lanes.  See \code{\link{betweenLaneNormalization}} for details.}

    \item{coerce}{\code{signature(from = "SeqExpressionSet", to = "CountDataSet")}: coercion to \code{DESeq} class \code{\linkS4class{CountDataSet}} for compatibility with downstream analysis.}
 }   
}

\author{Davide Risso <risso.davide@gmail.com>}

\seealso{
  \code{\linkS4class{eSet}}, \code{\link{newSeqExpressionSet}}, \code{\link{biasPlot}}, \code{\link{withinLaneNormalization}}, \code{\link{betweenLaneNormalization}}
}

\examples{

showMethods(class="SeqExpressionSet", where=getNamespace("EDASeq"))

exprs <- matrix(data=0,nrow=100,ncol=4)
for(i in 1:4) {
exprs[,i] <- rpois(100,lambda=50)
}
cond <- c(rep("A",2),rep("B",2))

counts <- newSeqExpressionSet(exprs,phenoData=AnnotatedDataFrame(data.frame(conditions=cond)))

head(exprs(counts))
boxplot(counts,col=as.numeric(pData(counts)[,1])+1)
}

\keyword{classes}
