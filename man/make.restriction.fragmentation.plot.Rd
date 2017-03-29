\name{make.restriction.fragmentation.plot}
\alias{make.restriction.fragmentation.plot}
\title{Plot the ratios of restriction fragmentation probes}
\description{Plot the ratios of restriction fragmentation probes, highlighting samples with low ratios }
\usage{
	make.restriction.fragmentation.plot(restr.data, low.ratio.samples = NULL)
}
\arguments{
  \item{restr.data}{A data-frame with raw counts of restriction fragmentation probes per sample. First row is mean count for probes A and B and second row is mean count for probes C and D.}
  \item{low.ratio.samples}{Samples previously identified as having low ratios}
}
\details{This function is called by \code{\link{restriction.fragmentation.qc}}.}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
