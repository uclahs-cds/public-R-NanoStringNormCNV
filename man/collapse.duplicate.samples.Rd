\name{collapse.duplicate.samples}
\alias{collapse.duplicate.samples}
\title{Collapse counts for duplicate samples}
\description{Collapse counts for duplicate samples. Only applicable if there are replicates}
\usage{
collapse.duplicate.samples(nano.df, annot)
}
\arguments{
  \item{nano.df}{A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{annot}{Annotations for the samples, must include the field Name}
}
\details{Collapse counts for duplicate samples. Only applicable if there are replicates. Collapsing is done with the mean count}
\value{A data-frame with collapsed counts for any duplicates}
\author{Cindy Yao}
\seealso{collapse.genes}
