\name{call.cnas.with.matched.normals}
\alias{call.cnas.with.matched.normals}
\title{Call CNAs using matched normal samples}
\description{Wrapper function that performs CNA calling on tumour-matched normal samples only}
\usage{
call.cnas.with.matched.normals(
  normalized.data, phenodata, per.chip = FALSE, call.method = 1,
  kd.values = c(0.85, 0.95), use.sex.info = TRUE
  )
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probe names must include "chrX" and "chrY" labels.}
  \item{phenodata}{Data-frame containing phenotype information. See output in \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use only reference sample(s) specific to each chip when calculating the tumour/normal ratios from which copy numbers are generated. Defaults to FALSE.}
  \item{call.method}{An integer specifying the CNA calling method passed on to \code{\link{call.copy.number.state}}:
  \cr1: round tumour/normal ratios to CNA calls using pre-defined thresholds derived from NanoString recommendations (Default)
  \cr2: use kernel density method to calculate copy number thresholds for CNA calling (see 'kd.values' argument)
  }
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values to be used to calculate copy number thresholds (when 'call.method' 2 is requested). Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain. For probes identified as belonging to male sex chromosomes (neutral CN = 1), deletions will be determined using the heterozygous threshold. Only required when kd.option = 3. Defaults to (0.85, 0.95), however, it is strongly recommended to customize these values for any given dataset.}
  \item{use.sex.info}{Whether to process sex chromosome probes separately. If TRUE (default), uses sample sex information in determining how to call copy number with these probes. It is not recommended to call CNAs on male sex chromosomes when there are few chrX/chrY probes per sample. See \code{\link{process.xy.probes}} for more information.}
}
\details{Get copy number states for each gene in each using corresponding matched normal samples. Be aware that tumour samples without a matched normal get dropped.}
\value{A list containing 2 items: 1) a gene X sample data frame of copy number states rounded according to the method specified 2) a gene X sample data frame of raw (not rounded) copy number states}
\author{Emilie Lalonde and Dorota Sendorek}
\seealso{call.cnas.with.pooled.normals}
