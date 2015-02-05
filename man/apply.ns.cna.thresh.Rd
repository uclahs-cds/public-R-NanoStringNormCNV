\name{apply.ns.cna.thresh}
\alias{apply.ns.cna.thresh}
\title{Determine the copy number state from counts}
\description{Determine the copy number state from counts using predefined thresholds as defined by NanoString}
\usage{
apply.ns.cna.thresh(tmr2ref)
}
\arguments{
  \item{tmr2ref}{The data-frame with tumour-normal signals to be converted to CNA states}
}
\details{Determine the copy number state for each gene in each sample based on predefined log2ratio thresholds.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
