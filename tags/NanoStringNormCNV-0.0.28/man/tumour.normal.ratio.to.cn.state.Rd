\name{tumour.normal.ratio.to.cn.state}
\alias{tumour.normal.ratio.to.cn.state}
\title{Convert tumour-normal ratio of CNA segments into discrete CNA calls}
\description{Will convert tumour-normal ratio CNA calls into CNA states for a set of probes}
\usage{tumour.normal.ratio.to.cn.state(ratio.data, thresholds)}
\arguments{
	\item{ratio.data}{A data frame of tumour-normal ratios per CNA segment/probe. Column must be labelled ratio}
	\item{thresholds}{A 4-item vector specifying homoz loss, het loss, single-copy gain and multi-level amplification, respectively. Order is important!)}
}
\details{Will convert tumour-normal ratio CNA calls into CNA states for a set of probes.}
\value{A vector of copy number states.}
\seealso{ \code{\link{get.tumour.normal.ratio}} }
\author{Emilie Lalonde}
