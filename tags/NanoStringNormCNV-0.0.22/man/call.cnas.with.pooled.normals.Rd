\name{call.cnas.with.pooled.normals}
\alias{call.cnas.with.pooled.normals}
\title{Call CNAs using a pooled reference}
\description{Wrapper function that performs CNA calling on tumour samples against a reference comprised of pooled normal samples.}
\usage{
call.cnas.with.pooled.normals(
  normalized.data, phenodata, per.chip = FALSE, call.method = 0, kd.values = NULL
  )
}
\arguments{
  \item{normalized.data}{The NanoString normalized counts for both test and control samples}
  \item{phenodata}{Data-frame containing phenotype information. See output in \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE.}
  \item{call.method}{An integer specifying the CNA calling method passed on to \code{\link{call.copy.number.state}}:
  0) round counts to CNA calls using NS pre-defined thresholds
  1) round counts to CNA calls using thresholds calculated from min/max values of normalized data
  2) use kernel density method to calculate CNA thresholds using default kernel density values (0.85, 0.95)
  3) use kernel density method to calculate CNA thresholds with user-provided kernel density values (see kd.values argument)
  Defaults to 0
  }
  \item{kd.values}{A 2- or 4-element vector containing kernel density values to use in defining thresholds for homozygous deletions, heterozygous deletion, single copy gain, multi-copy gain. Only required when kd.option = 3. Defaults to NULL}
}
\details{Get copy number states for each gene in each using corresponding matched normal samples. Preferable when most tumour samples do not have a matched normal. Pooled reference is created by taking the mean count of all normal samples.}
\value{A list containing 4 items (each being a gene X sample matrix):
1) copy number states rounded according to the method specified
2) raw (not rounded) copy number states
3) rounded copy number states for normal samples only
4) raw copy number states for normal samples only, with final column containing average ref CN}
\author{Emilie Lalonde and Dorota Sendorek}
