\name{tumour.normal.ratio.to.cn.state}
\alias{tumour.normal.ratio.to.cn.state}
\title{Convert tumour-normal ratios into discrete copy number (CN) states}
\description{Will convert tumour-normal count ratios into copy number states for a set of probes}
\usage{tumour.normal.ratio.to.cn.state(ratios, thresholds)}
\arguments{
	\item{ratios}{A vector of tumour/normal count ratios.}
	\item{thresholds}{A 4-item vector specifying homozygous loss, heterozygous loss, single-copy gain and multi-copy gain, respectively. Order is important!)}
}
\details{Converts tumour/normal ratios into CN states for a set of probes. Copy number states differ from actual copy numbers by representing the copy number deviation from an expected neutral copy number (i.e. 2): -2, -1, 0, 1, 2}
\value{A vector of copy number states.}
\seealso{ \code{\link{get.tumour.normal.ratio}}, \code{\link{get.sample.specific.cna.thresholds}} }
\author{Emilie Lalonde}
