\name{getGeneLengthAndGCContent}
\alias{getGeneLengthAndGCContent}

\title{Get gene length and GC-content}

\description{
Automatically retrieves gene length and GC-content information from Biomart or org.db packages.
}

\usage{
getGeneLengthAndGCContent(id, org, mode=c("biomart", "org.db"))
}

\arguments{
  \item{id}{
    Character vector of one or more ENSEMBL or ENTREZ gene IDs.
  }
  \item{org}{
    Organism three letter code, e.g. 'hsa' for 'Homo sapiens'. 
    See also: http://www.genome.jp/kegg/catalog/org_list.html;
    In org.db mode, this can be also a specific genome assembly, 
    e.g. 'hg38' or 'sacCer3'.
  }
  \item{mode}{
    Mode to retrieve the information. Defaults to 'biomart'. See Details.
  }
}

\details{
    The 'biomart' mode is based on functionality from the biomaRt 
    packgage and retrieves the required information from the BioMart database.
    This is available for all ENSEMBL organisms and is typically most current, 
    but can be time-consuming when querying several thousand genes at a time.

    The 'org.db' mode uses organism-based annotation packages from
    Bioconductor. This is much faster than the 'biomart' mode, but is only
    available for selected model organism currently supported by  
    BioC annotation functionality. 

    Results for the same gene ID(s) can differ between both modes as they
    are based on different sources for the underlying genome assembly. While
    the 'biomart' mode uses the latest ENSEMBL version, the 'org.db' mode uses
    BioC annotation packages typically built from UCSC. 
}

\value{
A numeric matrix with two columns: gene length and GC-content.
}

\seealso{
    \code{\link{getSequence}} to retrieve a genomic sequence from BioMart,
    \code{\link{genes}} to extract genomic coordinates from a TxDb object,
    \code{\link{getSeq}} to extract genomic sequences from a BSgenome object,
    \code{\link{alphabetFrequency}} to calculate nucleotide frequencies.
}

\author{
Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>
}

\examples{
getGeneLengthAndGCContent("ENSG00000012048", "hsa")
}
