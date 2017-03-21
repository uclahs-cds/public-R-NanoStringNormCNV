\name{call.cnas.with.pooled.normals}
\alias{call.cnas.with.pooled.normals}
\title{Call CNAs using a pooled normal reference}
\description{Wrapper function that performs CNA calling on tumour and normal samples against a reference comprised of pooled normal samples.}
\usage{
call.cnas.with.pooled.normals(
  normalized.data, phenodata, per.chip = FALSE, call.method = 0, 
  kd.values = NULL, use.sex.info = TRUE
  )
}
\arguments{
  \item{normalized.data}{A data-frame of normalized NanoString counts for both test and control samples. A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{Data-frame containing phenotype information. See output in \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE.}
  \item{call.method}{An integer specifying the CNA calling method passed on to \code{\link{call.copy.number.state}}:
  \cr0: round tumour/normal ratios to CNA calls using pre-defined thresholds derived from NanoString recommendations (Default)
  \cr1: round counts to CNA calls using thresholds calculated from min/max values of normalized data
  \cr2: use kernel density method to calculate CNA thresholds using default kernel density values (0.85, 0.95)
  \cr3: use kernel density method to calculate CNA thresholds with user-provided kernel density values (see kd.values argument)
  }
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values to be used to calculate copy number thresholds. Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain. For probes identified as belonging to male sex chromosomes (neutral CN = 1), deletions will be determined using the heterozygous threshold. Only required when kd.option = 3. Defaults to NULL}
  \item{use.sex.info}{Whether to process sex chromosome probes separately. If TRUE (default), uses sample sex information in determining how to call copy number with these probes. See \code{\link{process.xy.probes}} for more information.}
}
\details{Get copy number states for each gene in each using a reference of pooled normal samples. Preferable when most tumour samples do not have a matched normal. Pooled reference is created by taking the mean count of all normal samples.}
\value{A list containing 4 items (each being a gene X sample matrix):
\cr1) copy number states rounded according to the method specified
\cr2) raw (not rounded) copy number states
\cr3) rounded copy number states for normal samples only
\cr4) raw copy number states for normal samples only, with final column containing average ref CN}
\author{Emilie Lalonde and Dorota Sendorek}
\seealso{call.cnas.with.matched.normals}
