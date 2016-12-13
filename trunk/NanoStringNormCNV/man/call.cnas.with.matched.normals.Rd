\name{call.cnas.with.matched.normals}
\alias{call.cnas.with.matched.normals}
\title{Call CNAs using matched normal samples}
\description{Wrapper function that performs CNA calling on tumour-matched normal samples only}
\usage{
call.cnas.with.matched.normals(
  normalized.data, phenodata, per.chip = FALSE, call.method = 1,
  kd.values = NULL, use.sex.info = TRUE
  )
}
\arguments{
  \item{normalized.data}{A data-frame of normalized NanoString counts for both test and control samples. A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{Data-frame containing phenotype information. See output in \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE.}
  \item{call.method}{An integer specifying the CNA calling method passed on to \code{\link{call.copy.number.state}}:
  1) round counts to CNA calls using NS pre-defined thresholds
  2) use kernel density method to calculate CNA thresholds using default kernel density values (0.85, 0.95)
  3) use kernel density method to calculate CNA thresholds with user-provided kernel density values (see kd.values argument)
  Defaults to 1
  }
  \item{kd.values}{A 2- or 4-element vector containing kernel density values to use in defining thresholds for homozygous deletions, heterozygous deletion, single copy gain, multi-copy gain. Only required when kd.option = 3. Defaults to NULL}
  \item{use.sex.info}{Whether to process sex chromosome probes separately. If TRUE (default), uses sample sex information in determining how to call copy number with these probes. It is not recommended to call CNAs on male sex chromosomes when there are few chrX/chrY probes per sample. See \code{\link{process.xy.probes}} for more information.}
}
\details{Get copy number states for each gene in each using a reference of pooled normal samples. Be aware that tumour samples without a matched normal get dropped.}
\value{A list containing 2 items: 1) a gene X sample data frame of copy number states rounded according to the method specified 2) a gene X sample data frame of raw (not rounded) copy number states}
\author{Emilie Lalonde and Dorota Sendorek}
