\name{make.restr.digest.plot}
\alias{make.restr.digest.plot}
\title{Plot the ratio of restriction digestion probes}
\description{Plot the ratio of restriction digestion probes and highlight those with low-ratios}
\usage{
make.restr.digest.plot(restr.data, low.ratio.samples = NULL)
}
\arguments{
  \item{restr.data}{A data-frame with the raw counts for the restriction digestion probes}
  \item{low.ratio.samples}{Samples previously identified as having low-ratios}
}
\details{Should be called from within restr.digest.norm}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
