\name{make.invariant.probe.plot}
\alias{make.invariant.probe.plot}
\title{Make diagnostic plots for invariant probes}
\description{Make diagnostic plots for invariant probes}
\usage{
make.invariant.probe.plot(inv.probe.counts, tissue.type = NULL)
}
\arguments{
  \item{inv.probe.counts}{Counts for invariant probes}
  \item{tissue.type}{Data frame containing 'SampleID' (first column) and 'Type' -- 'Tumour' or 'Reference' (second column). Defaults to NULL}
}
\details{Creates a scatterplot showing the counts for each invariant probe, and a barplot showing the number of probes with low coverage per sample}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
