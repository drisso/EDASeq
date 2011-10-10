\name{EDASeq-package}
\alias{EDASeq-package}
\alias{EDASeq}
\docType{package}
\title{Exploratory Data Analysis and Normalization for RNA-Seq data}
\description{Numerical summaries and graphical representations of some key features of the data along with implementations of both within-lane normalization methods for GC content bias and between-lane normalization methods to adjust for sequencing depth and possibly other differences in distribution.}
\details{
The \code{\linkS4class{SeqExpressionSet}} class is used to store gene-level counts along with sample information. It extends the virtual class \code{\linkS4class{eSet}}. See the help page of the class for details. 

"Read-level" information is managed via the \code{\linkS4class{FastqFileList}} and \code{\linkS4class{BamFileList}} classes of \code{\link{Rsamtools}}.

Most used graphic tools for the \code{\linkS4class{FastqFileList}} and \code{\linkS4class{BamFileList}} objects are: 'barplot', 'plotQuality', 'plotNtFrequency'. For \code{\linkS4class{SeqExpressionSet}} objects are: 'biasPlot', 'meanVarPlot', 'MDPlot'.

To perform gene-level normalization use the functions 'withinLaneNormalization' and 'betweenLaneNormalization'.

An 'As' method exists to coerce \code{\linkS4class{SeqExpressionSet}} objects to \code{\linkS4class{CountDataSet}} objects (\code{DESeq} package).

See the package vignette for a typical Exploratory Data Analysis example.
}
\author{
Davide Risso and Sandrine Dudoit.
Maintainer: Davide Risso <risso.davide@gmail.com>
}
\references{
J. H. Bullard, E. A. Purdom, K. D. Hansen and S. Dudoit (2010). Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics Vol. 11, Article 94.

D. Risso, K. Schwartz, G. Sherlock and S. Dudoit (2011). GC-Content Normalization for RNA-Seq Data. Technical Report No. 291, Division of Biostatistics, University of California, Berkeley, Berkeley, CA.
}
