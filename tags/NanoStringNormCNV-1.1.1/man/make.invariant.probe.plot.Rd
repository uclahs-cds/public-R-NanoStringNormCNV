\name{make.invariant.probe.plot}
\alias{make.invariant.probe.plot}
\title{Make diagnostic plots for invariant probes}
\description{Make diagnostic plots for invariant probes}
\usage{
	make.invariant.probe.plot(inv.probe.counts, tissue.type = NULL)
}
\arguments{
  \item{inv.probe.counts}{A gene by sample data-frame of counts from invariant probes}
  \item{tissue.type}{A data-frame of sample annotation with columns 'SampleID' and 'Type' (specifying whether 'Tumour' or 'Reference'). Recommended to use output of \code{\link{load.phenodata}}. Defaults to NULL}
}
\details{
	Creates a scatterplot showing the counts for each invariant probe and a barplot showing the number of probes with low coverage per sample (if any)
}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
