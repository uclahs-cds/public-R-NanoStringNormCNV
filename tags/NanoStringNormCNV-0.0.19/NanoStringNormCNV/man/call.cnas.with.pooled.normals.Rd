\name{call.cnas.with.pooled.normals}
\alias{call.cnas.with.pooled.normals}
\title{Call CNAs using a pooled reference}
\description{Wrapper function that performs CNA calling on tumour samples against a reference comprised of pooled normal samples.}
\usage{
call.cnas.with.pooled.normals(normalized.data, phenodata, per.chip = FALSE, kd.option = 0, kd.values = NULL)
}
\arguments{
  \item{normalized.data}{The NanoString normalized counts for both test and control samples}
  \item{phenodata}{Data-frame containing phenotype information. See output in \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE.}
  \item{kd.option}{An integer specifying the method passed on to 'call.copy.number.state': 0) segment using NS pre-defined thresholds (Default)
  1) define thresholds based on values observed in reference samples
  2) create thresholds with kernel density approach using default KD values
  3) create thresholds with kernel density approach using user-provided KD values }
  \item{kd.values}{KD values to use in defining thresholds. Required when kd.option = 3. Defaults to NULL}
}
\details{Get copy number states for each gene in each using corresponding matched normal samples. Preferable when most tumour samples do not have a matched normal. Pooled reference is created by taking the mean count.}
\value{A list containing 4 items (each being a gene X sample matrix):
1) copy number states rounded according to the method specified
2) raw (not rounded) copy number states
3) rounded copy number states for normal samples only
4) raw copy number states for normal samples only, with final column containing average ref CN}
\author{Emilie Lalonde and Dorota Sendorek}
