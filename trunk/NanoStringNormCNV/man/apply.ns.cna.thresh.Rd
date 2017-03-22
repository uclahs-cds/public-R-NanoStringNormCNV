\name{apply.ns.cna.thresh}
\alias{apply.ns.cna.thresh}
\title{Determine the copy number (CN) values from pre-defined count thresholds}
\description{Determine the copy number values from count ratios using pre-defined thresholds}
\usage{
apply.ns.cna.thresh(ratio.data, cna.thresh = c(0.4, 1.5, 2.5, 3.5))
}
\arguments{
  \item{ratio.data}{A gene by sample data-frame of tumour/normal ratios. See output of \code{\link{get.tumour.normal.ratio}}}
  \item{cna.thresh}{A 4-element vector with the tumour-normal thresholds to call homozygous deletions, heterozygous deletions, single copy gains, and multi-copy gains. Default values (0.4, 1.5, 2.5, 3.5) are derived from NanoString recommendations.}
}
\details{Determine the actual copy number value (0, 1, 2, 3, 4) for each probe using pre-defined thresholds for tumour/normal count ratios.}
\value{A gene by sample data-frame with copy number values}
\author{Cindy Yao and Emilie Lalonde}
