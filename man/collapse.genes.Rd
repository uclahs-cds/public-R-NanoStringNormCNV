\name{collapse.genes}
\alias{collapse.genes}
\title{Collapse counts for probes in the same gene}
\description{Collapse counts for probes in the same gene. Only applicable if there are multiple probes per gene}
\usage{
collapse.genes(nano.df)
}
\arguments{
  \item{nano.df}{The NanoString data-frame with normalized counts}
}
\details{Collapse counts for probes within the same gene. Collapsing is done with the mean count}
\value{A data-frame with collapsed counts}
\author{Cindy Yao}
\seealso{collapse.duplicate.samples}
