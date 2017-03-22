\name{collapse.genes}
\alias{collapse.genes}
\title{Collapse counts for probes in the same gene}
\description{Collapse counts of all probes belonging to the same gene}
\usage{
	collapse.genes(normalized.data)
}
\arguments{
  \item{normalized.data}{A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Probes to be collapsed together are identified by matching 'Accession' values.}
}
\details{Collapse counts for probes belonging to the same gene (only performed on 'Endogenous' or 'Housekeeping' probes). Collapsing is done by taking the mean count.}
\value{A gene by sample data-frame of collapsed counts, where first three columns contain 'CodeClass', 'Name', and 'Accession' information}
\author{Cindy Yao}
\seealso{\code{\link{collapse.replicates}}}
