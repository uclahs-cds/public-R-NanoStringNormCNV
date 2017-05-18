\name{apply.kd.cna.thresh}
\alias{apply.kd.cna.thresh}
\title{Determine the copy number (CN) values from ratios using a kernel density approach}
\description{Determine the copy number values from count ratios by first calculating thresholds using a kernel density approach}
\usage{
apply.kd.cna.thresh(ratio.data, kd.values, neutral.cn = 2)
}
\arguments{
  \item{ratio.data}{A gene by sample data-frame of tumour:normal ratios. See output of \code{\link{get.tumour.normal.ratio}}}
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values to be used to calculate copy number thresholds. Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain.}
  \item{neutral.cn}{The expected neutral copy number; i.e. 2 for autosomes/female chrX (default), 1 for male sex chromosomes. This value is used to calculate the observed copy number values by summing it with generated copy number states: -2, -1, 0, 1, 2 (see output of \code{\link{tumour.normal.ratio.to.cn.state}}).}
}
\details{Determine the copy number values (0, 1, 2, 3, 4) of each probe from raw tumour:normal count ratios using thresholds derived from requested kernel densities ('kd.values').}
\value{A gene by sample data-frame of copy number values}
\author{Cindy Yao and Emilie Lalonde}
\seealso{\code{\link{apply.ns.cna.thresh}}}
