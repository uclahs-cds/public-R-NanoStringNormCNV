\name{tumour.normal.ratio.to.cn.state}
\alias{tumour.normal.ratio.to.cn.state}
\title{Convert tumour:normal count ratios into discrete copy number (CN) states}
\description{Convert tumour:normal count ratios into copy number states for a set of probes}
\usage{tumour.normal.ratio.to.cn.state(ratios, thresholds)}
\arguments{
	\item{ratios}{A vector of tumour:normal count ratios. See output of \code{\link{get.tumour.normal.ratio}}}
	\item{thresholds}{A 4-item vector specifying homozygous loss, heterozygous loss, single-copy gain and multi-copy gain, respectively. Order is important! See output of \code{\link{get.cna.thresholds}}}
}
\details{Copy number states differ from actual copy number values by representing the copy number deviation (i.e. -2, -1, 0, 1, 2+) from an expected neutral copy number (i.e. 2) and not the copy number itself (i.e. 0, 1, 2, 3, 4+).}
\value{A vector of copy number states}
\author{Emilie Lalonde}
