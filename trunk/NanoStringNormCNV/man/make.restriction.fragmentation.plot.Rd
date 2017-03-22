\name{make.restriction.fragmentation.plot}
\alias{make.restriction.fragmentation.plot}
\title{Plot the ratios of restriction fragmentation probes}
\description{Plot the ratios of restriction fragmentation probes, highlighting samples with low ratios}
\usage{
	make.restriction.fragmentation.plot(restr.data, low.ratio.samples = NULL)
}
\arguments{
  \item{restr.data}{A gene by sample data-frame of raw NanoString counts of restriction fragmentation probes. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{low.ratio.samples}{Samples previously identified as having low ratios}
}
\details{This function is called by \code{\link{restriction.fragmentation.qc}}.}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
