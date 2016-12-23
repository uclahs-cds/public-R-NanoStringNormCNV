\name{apply.kd.cna.thresh}
\alias{apply.kd.cna.thresh}
\title{Determine the copy number state from counts}
\description{Determine the copy number state from counts using a kernel density approach}
\usage{
apply.kd.cna.thresh(ratio.data, kd.values)
}
\arguments{
  \item{ratio.data}{A gene by sample data-frame of tumour/normal ratios to be converted to CNA states}
  \item{kd.values}{A 2- or 4-element vector containing the KD thresholds to be used for (hom, het) deletions and (het, hom) amplifications}
}
\details{Determine the copy number state for each gene in each sample based on kernel densities using the provided proportions.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
