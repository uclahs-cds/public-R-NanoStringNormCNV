\name{make.invariant.probe.plot}
\alias{make.invariant.probe.plot}
\title{Make diagnostic plots for invariant probes}
\description{Make diagnostic plots for invariant probes}
\usage{
make.invariant.probe.plot(inv.probe.counts, fname.stem)
}
\arguments{
  \item{inv.probe.counts}{Counts for invariant probes}
  \item{fname.stem}{Name of type of invariant probes (invariant_probes or housekeeping_probes or similar is suggested)}
}
\details{Creates a scatterplot showing the counts for each invariant probe, and a barplot showing the number of probes with low coverage per sample}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
