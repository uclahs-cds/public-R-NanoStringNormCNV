\name{make.invariant.probe.plot}
\alias{make.invariant.probe.plot}
\title{Make diagnostic plots for invariant probes}
\description{Make diagnostic plots for invariant probes}
\usage{
	make.invariant.probe.plot(inv.probe.counts, tissue.type = NULL)
}
\arguments{
  \item{inv.probe.counts}{A gene by sample data-frame of counts from invariant probes}
  \item{tissue.type}{A data-frame with fields 'SampleID' and 'Type' (specifying whether 'Tumour' or 'Reference'). Defaults to NULL}
}
\details{
	Creates a scatterplot showing the counts for each invariant probe and a barplot showing the number of probes with low coverage per sample (if applicable)
}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
