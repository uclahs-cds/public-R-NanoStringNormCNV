\name{collapse.replicates}
\alias{collapse.replicates}
\title{Collapse probe counts of sample replicates}
\description{Collapse probe counts of sample replicates}
\usage{
	collapse.replicates(normalized.data, phenodata)
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{A data-frame of sample annotation which includes the fields 'SampleID' and 'Name'. Sample replicates are identified by their matching 'Name' values. Recommended to use output of \code{\link{load.phenodata}}}
}
\details{Collapsing is done by taking the median count.}
\value{A gene by sample data-frame of collapsed counts, where first three column names are 'CodeClass', 'Name', and 'Accession'.}
\author{Cindy Yao}
\seealso{\code{\link{collapse.genes}}}
