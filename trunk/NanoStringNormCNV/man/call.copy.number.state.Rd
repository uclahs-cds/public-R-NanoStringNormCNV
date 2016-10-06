\name{call.copy.number.state}
\alias{call.copy.number.state}
\title{Determine the copy number state from counts}
\description{Subtract the reference signal from the test sample}
\usage{
call.copy.number.state(
  input, reference, per.chip = FALSE, chip.info = NULL,
  thresh.method = 'round', multi.factor = 2, kd.vals = c(0.85, 0.95),
  adjust = FALSE, cna.thresh = c(0.4, 1.5, 2.5, 3.5)
  )
}
\arguments{
  \item{input}{The NanoString normalized counts for both test and control samples}
  \item{reference}{The name(s) of reference sample(s) to be used. Must match the column header of input}
  \item{per.chip}{Use the reference specific to each chip only. Defaults to FALSE.}
  \item{chip.info}{Annotations for the chips, must include the fields Chip and SampleID. Defaults to NULL}
  \item{thresh.method}{How to assign CNA states- default is the NS recommended method using their thresholds ('round'), see cna.thresh for control. Can also specify 'KD' for kernel density, see kd.vals for control. Otherwise specify 'none'}
  \item{multi.factor}{A number to multiply the test/ref ratio by. NanoString recommends 2 (default), otherwise use 1 for no amplification/boosting}
  \item{kd.vals}{A 2- or 4-element vector containing the KD thresholds to be used for (hom, het) deletions and (het, hom) amplifications. Defaults to 0.85, 0.95}
  \item{adjust}{Whether or not to center the data at 2 prior to calling CNAs, default is FALSE}
  \item{cna.thresh}{A 4-element vector with the tumour-normal thresholds to call homoz deletions, het deletion, single copy gain, multi-copy gain. The defaults are the recommended values by NanoString: 0.4, 1.5, 2.5, 3.5}
}
\details{Determine the copy number state for each gene in each sample based on counts in the test and reference set.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
