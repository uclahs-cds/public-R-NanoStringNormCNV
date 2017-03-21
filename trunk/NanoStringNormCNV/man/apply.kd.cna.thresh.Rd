\name{apply.kd.cna.thresh}
\alias{apply.kd.cna.thresh}
\title{Determine the copy number (CN) states by calculating count thresholds}
\description{Determine the copy number state from counts using a kernel density approach}
\usage{
apply.kd.cna.thresh(ratio.data, kd.values)
}
\arguments{
  \item{ratio.data}{A gene by sample data-frame of tumour/normal ratios to be converted to CNA states. See output of \code{\link{get.tumour.normal.ratio}}}
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values to be used to calculate copy number thresholds. Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain.}
}
\details{Determine the copy number state (-2, -1, 0, 1, 2+) of each probe using thresholds derived from the kernel densities of requested proportions ('kd.values'). Resulting states will have to be summed with the expected neutral copy number (i.e. 2) to get observed copy number values.}
\value{A gene by sample data-frame of copy number states}
\author{Cindy Yao and Emilie Lalonde}
