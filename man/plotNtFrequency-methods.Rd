\name{plotNtFrequency-methods}
\docType{methods}
\alias{plotNtFrequency}
\alias{plotNtFrequency-methods}
\alias{plotNtFrequency,BamFile-method}
\alias{plotNtFrequency,ShortRead-method}
\title{ Methods for Function \code{plotNtFrequency} in Package \pkg{EDASeq} }
\description{
Plots the nucleotide frequencies per position.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "ShortRead")}}{
}

\item{\code{signature(x = "BamFile")}}{
}

}

It plots the nucleotide frequencies per position, averaging all the reads in \code{x}.

}
\keyword{methods}

\examples{
if(require(leeBamViews)) {
  files <- list.files(file.path(system.file(package = "leeBamViews"), "bam"),
                    pattern = "bam$", full.names = TRUE)
  names(files) <- gsub("\\.bam", "", basename(files))

  bfs <- BamFileList(files)

  plotNtFrequency(bfs[[1]])
}
}

