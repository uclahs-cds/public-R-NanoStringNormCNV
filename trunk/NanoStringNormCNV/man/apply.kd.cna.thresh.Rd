\name{apply.kd.cna.thresh}
\alias{apply.kd.cna.thresh}
\title{Determine the copy number state from counts}
\description{Determine the copy number state from counts using a kernel density approach}
\usage{
apply.kd.cna.thresh(tmr2ref, kd.thresh, sex.info, sex.probes = NULL)
}
\arguments{
  \item{tmr2ref}{The data-frame with tumour-normal signals to be converted to CNA states}
  \item{kd.thresh}{A 2- or 4-element vector containing the KD thresholds to be used for (hom, het) deletions and (het, hom) amplifications}
  \item{sex.info}{Data frame containing columns 'SampleID' and 'sex' for each sample in input}
  \item{sex.probes}{A vector of names of sex chromosome probes found in dataset. Required if any samples in tmr2ref are male. Defaults to NULL}
}
\details{Determine the copy number state for each gene in each sample based on kernel densities using the provided proportions.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
