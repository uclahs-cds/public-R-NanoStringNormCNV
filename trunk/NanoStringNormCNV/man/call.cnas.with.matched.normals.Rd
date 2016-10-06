\name{call.cnas.with.matched.normals}
\alias{call.cnas.with.matched.normals}
\title{Call CNAs using matched normal samples}
\description{Wrapper function that performs CNA calling on tumour-matched normal samples only}
\usage{
call.cnas.with.matched.normals(
  normalized.data, phenodata, per.chip = FALSE, kd.option = 0
  )
}
\arguments{
  \item{normalized.data}{The NanoString normalized counts for both test and control samples}
  \item{phenodata}{Data-frame containing phenotype information. See output in \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE.}
  \item{kd.option}{An integer specifying the method passed on to \code{\link{call.copy.number.state}}:
  0|1) segment using NS pre-defined thresholds
  2|3) segment using kernel densitys
  Defaults to 0
  }
}
\details{Get copy number states for each gene in each using a reference of pooled normal samples. Be aware that tumour samples without a matched normal get dropped.}
\value{A list containing 2 items: 1) a gene X sample matrix of copy number states rounded according to the method specified 2) a gene X sample matrix of raw (not rounded) copy number states}
\author{Emilie Lalonde and Dorota Sendorek}
