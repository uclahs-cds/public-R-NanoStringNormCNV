\name{apply.ns.cna.thresh}
\alias{apply.ns.cna.thresh}
\title{Determine the copy number state from counts}
\description{Determine the copy number state from count ratios using pre-defined thresholds as suggested by NanoString}
\usage{
apply.ns.cna.thresh(ratio.data, cna.thresh = c(0.4, 1.5, 2.5, 3.5))
}
\arguments{
  \item{ratio.data}{A gene by sample data-frame of tumour/normal ratios to be converted to CNA states}
  \item{cna.thresh}{A 4-element vector with the tumour-normal thresholds to call homoz deletions, het deletion, single copy gain, multi-copy gain. The defaults are the recommended values by NanoString}
}
\details{Determine the copy number state for each gene in each sample based on pre-defined log2ratio thresholds.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
