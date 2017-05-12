\name{call.copy.number.values}
\alias{call.copy.number.values}
\title{Determine the copy number from counts}
\description{Uses tumour/normal count ratios to determine segment copy number (CN)}
\usage{
call.copy.number.values(
  normalized.data, reference, per.chip = FALSE, chip.info = NULL,
  thresh.method = 'round', cna.thresh = c(0.4, 1.5, 2.5, 3.5), kd.values = c(0.85, 0.95),
  multi.factor = 2, adjust = FALSE 
  )
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probes are identified by "chrX" and "chrY" labels in their 'Name' and/or 'Accession' values.}
  \item{reference}{Name(s) of reference sample(s) to be used. Must match the sample ID column names in 'normalized.data'}
  \item{per.chip}{Whether to use the reference sample(s) specific to a given chip only. Defaults to FALSE.}
  \item{chip.info}{A data-frame of chip information for each sample in 'normalized.data', including columns 'Cartridge' and 'SampleID', or NULL (default). Required if 'per.chip' = TRUE. Recommended to use output of \code{\link{load.phenodata}}}
  \item{thresh.method}{How to determine copy numbers. Default is 'round', where calling is performed using pre-defined thresholds (see 'cna.thresh' for control). Alternatively, can use 'KD', where dataset-specific thresholds are first obtained using a kernel density approach (see 'kd.values' for control). To get raw tumour/normal ratios, use 'none'.}
  \item{cna.thresh}{A 4-element vector of tumour/normal thresholds with which to call homozygous deletions, heterozygous deletions, single copy gains, multi-copy gains. Default values are derived from NanoString recommendations: 0.4, 1.5, 2.5, 3.5.}
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values (between 0 and 1) to be used to calculate copy number thresholds (when 'call.method' 2 is requested). Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain. For probes identified as belonging to male sex chromosomes (neutral CN = 1), deletions will be determined using the heterozygous threshold. Default values are: 0.85, 0.95 (alternatively, this can be written as: 0.95, 0.85, 0.85, 0.95). However, it is **strongly recommended** to tune this parameter to any given dataset.}
  \item{multi.factor}{A number by which to multiply the raw tumour/normal ratios prior to CN calling. Also, the expected neutral CN used to calculate observed CN values if 'thresh.method' = 'KD' (in this scenario, CN states (-2, -1, 0, 1, 2) are generated and used to determine actual CN by summing neutral CN and CN state together. See \code{\link{apply.kd.cna.thresh}}). NanoString recommends a value of 2 (default), otherwise use 1 for no amplification/boosting (i.e. with chrX and chrY probes in male samples).}
  \item{adjust}{Whether or not to center the data around the 'multi.factor' value (i.e. 2) prior to calling CNAs. Default is FALSE}
}
\details{Determine the copy number values for each probe using the count differences between the test and reference sets.}
\references{
  See NanoString website for PDFs on analysis guidelines:
    \url{https://www.nanostring.com/support/product-support/support-documentation}

  The NanoString assay is described in the paper:
    Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\value{A gene by sample data-frame of copy number values}
\author{Cindy Yao and Emilie Lalonde}
