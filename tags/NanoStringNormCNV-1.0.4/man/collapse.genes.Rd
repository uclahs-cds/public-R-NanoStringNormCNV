\name{collapse.genes}
\alias{collapse.genes}
\title{Collapse probe counts per gene}
\description{Collapse counts of all probes belonging to the same gene}
\usage{
	collapse.genes(normalized.data)
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Probes to be collapsed together are identified by their matching 'Accession' values.}
}
\details{Collapse counts for probes belonging to the same gene (only performed on 'Endogenous' or 'Housekeeping' probes). Collapsing is done by taking the mean count.}
\value{A gene by sample data-frame of collapsed counts, where first three column names are 'CodeClass', 'Name', and 'Accession'.}
\author{Cindy Yao}
