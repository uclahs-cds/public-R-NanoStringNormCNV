\name{collapse.replicates}
\alias{collapse.replicates}
\title{Collapse counts for sample replicates}
\description{Collapse counts of sample replicates}
\usage{
	collapse.replicates(normalized.data, phenodata)
}
\arguments{
  \item{normalized.data}{A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{A data-frame of sample annotation which includes the fields 'SampleID' and 'Name'. Sample replicates are identified by their matching 'Name' values. See output in \code{\link{load.phenodata}}}
}
\details{Collapse counts of sample replicates by taking the median count.}
\value{A gene by sample data-frame of collapsed counts, where first three columns contain 'CodeClass', 'Name', and 'Accession' information}
\author{Cindy Yao}
\seealso{collapse.genes}
