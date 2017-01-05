\name{call.copy.number.state}
\alias{call.copy.number.state}
\title{Determine the copy number from counts}
\description{Uses tumour/normal count ratios to determine segment copy number}
\usage{
call.copy.number.state(
  normalized.data, reference, per.chip = FALSE, chip.info = NULL,
  thresh.method = 'round', cna.thresh = c(0.4, 1.5, 2.5, 3.5), kd.values = c(0.85, 0.95),
  multi.factor = 2, adjust = FALSE 
  )
}
\arguments{
  \item{normalized.data}{The NanoString normalized counts for both test and control samples. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probe names must include "chrX" and "chrY" labels.}
  \item{reference}{The name(s) of reference sample(s) to be used. Must match the column names of 'normalized.data'}
  \item{per.chip}{Whether to use the reference sample(s) specific to a given chip only. Defaults to FALSE.}
  \item{chip.info}{A data-frame of chip information for each sample in 'normalized.data' with columns 'Cartridge' and 'SampleID' or NULL (default). Required if 'per.chip' = TRUE. See output of \code{\link{load.phenodata}}}
  \item{thresh.method}{How to determine copy number. Default value is 'round', where calling is performed using pre-defined thresholds, see 'cna.thresh' for control. Can also use 'KD', where dataset-specific thresholds are first calculated using a kernel density approach, see 'kd.values' for control. Otherwise, specify 'none' to get raw tumour/normal ratios.}
  \item{cna.thresh}{A 4-element vector with the tumour-normal thresholds to call homozygous deletions, heterozygous deletions, single copy gains, multi-copy gains. Default values are derived from NanoString recommendations: (0.4, 1.5, 2.5, 3.5)}
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values to be used to calculate (homozygous, heterozygous) deletion and (heterozygous, homozygous) amplification thresholds. Default values are (0.85, 0.95). However, it is **strongly recommended** to fit this parameter to a given dataset.}
  \item{multi.factor}{A number to multiply the tumour/normal ratios by. Also, the expected neutral copy number used in calculating observed copy numbers when 'thresh.method' = 'KD' (in this scenario, copy number states (-2, -1, 0, 1, 2) are first generated and then used to determine actual copy number by summing neutral CN and CN state together. See \code{\link{apply.kd.cna.thresh}}). NanoString recommends 2 (default), otherwise use 1 for no amplification/boosting (i.e. with chrX and chrY probes in male samples).}
  \item{adjust}{Whether or not to center the data at neutral copy number (i.e. 2) prior to calling CNAs, default is FALSE}
}
\details{Determine the copy number for each gene in each sample based on count differences between the test and reference sets.}
\value{A gene by sample data-frame of copy number values}
\author{Cindy Yao and Emilie Lalonde}
