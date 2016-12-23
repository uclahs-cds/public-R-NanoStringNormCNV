\name{apply.kd.cna.thresh}
\alias{apply.kd.cna.thresh}
\title{Determine the copy number (CN) states by calculating count thresholds}
\description{Determine the copy number state from counts using a kernel density approach}
\usage{
apply.kd.cna.thresh(ratio.data, kd.values)
}
\arguments{
  \item{ratio.data}{A gene by sample data-frame of tumour/normal ratios to be converted to CNA states. See output of \code{\link{get.tumour.normal.ratio}}}
  \item{kd.values}{A 2- or 4-element vector containing the kernel density thresholds to be used for (homozygous, heterozygous) deletions and (heterozygous, homozygous) amplifications}
}
\details{Determine the copy number state (-2, -1, 0, 1, 2) for each gene in each sample based on kernel densities using the provided proportions.}
\value{A gene by sample data-frame of copy number states}
\author{Cindy Yao and Emilie Lalonde}
