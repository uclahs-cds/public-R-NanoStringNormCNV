\name{collapse.genes}
\alias{collapse.genes}
\title{Collapse counts for probes in the same gene}
\description{Collapse counts for probes in the same gene. Only applicable if there are multiple probes per gene}
\usage{
collapse.genes(normalized.data)
}
\arguments{
  \item{normalized.data}{A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
}
\details{Collapse counts for 'Endogenous' and 'Housekeeping' probes within the same gene. Collapsing is done with the mean count}
\value{A data-frame with collapsed counts}
\author{Cindy Yao}
\seealso{collapse.duplicate.samples}
