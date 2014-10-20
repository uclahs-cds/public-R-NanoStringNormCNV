\name{call.copy.number.state}
\alias{call.copy.number.state}
\title{Determine the copy number state from counts}
\description{Subtract the reference signal from the test sample}
\usage{
call.copy.number.state(
	input, reference, per.chip = FALSE, chip.info = NULL, to.round = TRUE, to.log = FALSE,
	multi.factor = 2
	)
}
\arguments{
  \item{input}{The NanoString normalized counts for both test and control samples}
  \item{reference}{The name(s) of reference sample(s) to be used. Must match the column header of input}
  \item{per.chip}{Use the reference specific to each chip only. Default to false.}
  \item{chip.info}{Annotations for the chips, must include the fields Chip and SampleID}
  \item{to.round}{Whether to round the copy number states to the nearest integer. Default to true.}
  \item{to.log}{Whether to log the normalized data prior to copy number state calculation.}
  \item{multi.factor}{A number to multiply the test/ref ratio. NanoString recommends 2 (default), otherwise use 1 for no amplification/boosting}
}
\details{Determine the copy number state for each gene in each sample based on counts in the test and reference set.}
\value{A data-frame with copy number states}
\author{Cindy Yao and Emilie Lalonde}
\seealso{aggregate.genes}
