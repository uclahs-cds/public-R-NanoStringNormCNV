\name{apply.ns.cna.thresh}
\alias{apply.ns.cna.thresh}
\title{Determine the copy number state from counts}
\description{Determine the copy number state from counts using predefined thresholds as defined by NanoString}
\usage{
apply.ns.cna.thresh(tmr2ref, thresh = c(0.4, 1.5, 2.5, 3.5))
}
\arguments{
  \item{tmr2ref}{The data-frame with tumour-normal signals to be converted to CNA states}
  \item{thresh}{A vector of length 4 with the tumour-normal thresholds to call homoz deletions, het deletion, single copy gain, multi-copy gain. The defaults are the recommended values by NS}
}
\details{Determine the copy number state for each gene in each sample based on predefined log2ratio thresholds.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
