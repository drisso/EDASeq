\name{plotQuality-methods}
\docType{methods}
\alias{plotQuality}
\alias{plotQuality-methods}
\alias{plotQuality,BamFileList-method}
\alias{plotQuality,BamFile-method}
\alias{plotQuality,FastqFileList-method}
\title{  Methods for Function \code{plotQuality} in Package \pkg{EDASeq} }
\description{
\code{plotQuality} produces a plot of the quality of the reads.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "BamFileList")}}{
It produces a plot that summarizes the per-base mean quality of the reads of each BAM file in \code{x}.
}

\item{\code{signature(x = "BamFile")}}{
It produces a boxplot of the per-base distribution of the quality scores of the reads in \code{x}.
}

\item{\code{signature(x = "FastqFileList")}}{
It produces a plot that summarizes the per-base mean quality of the reads of each FASTQ file in \code{x}. 
}

}}

\section{Details}{
Since FASTQ files can be very long, it can be very expensive to process a whole file. One way to avoid this, is to consider a subset of the file and then plot the quality of the subset. As long as one assumes that the subset is random, this is a good approximation. The function \code{\link{FastqSampler}} of \code{ShortRead} can be used for this. See its help page for an example.
}

\keyword{methods}



