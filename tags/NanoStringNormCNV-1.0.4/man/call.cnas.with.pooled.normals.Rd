\name{call.cnas.with.pooled.normals}
\alias{call.cnas.with.pooled.normals}
\title{Call CNAs using a pooled normal reference}
\description{Wrapper function that performs CNA calling on both tumour and normal samples against a reference comprised of pooled normal samples.}
\usage{
call.cnas.with.pooled.normals(
  normalized.data, phenodata, per.chip = FALSE, call.method = 1, 
  kd.values = c(0.85, 0.95), use.sex.info = TRUE
  )
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probes are identified by "chrX" and "chrY" labels in their 'Name' and/or 'Accession' values.}
  \item{phenodata}{A data-frame of sample annotation. Recommended to use output of \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use reference samples specific to each chip when calculating the tumour-normal ratios from which copy numbers are generated. Defaults to FALSE.}
  \item{call.method}{The CNA calling method passed on to \code{\link{call.copy.number.values}}. Options are:
  \cr1: round tumour/normal ratios to CNA calls using pre-defined thresholds derived from NanoString recommendations (Default)
  \cr2: use kernel density method to calculate copy number thresholds for CNA calling (see 'kd.values' argument)
  \cr3: round tumour/normal ratios to CNA calls using thresholds calculated from min/max values of normal sample counts (offset = (median max - median min) * 0.15)
  }
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values (between 0 and 1) to be used to calculate copy number thresholds (when 'call.method' 2 is requested). Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain. For probes identified as belonging to male sex chromosomes (neutral CN = 1), deletions will be determined using the heterozygous threshold. Default values are: 0.85, 0.95 (alternatively, this could be written as: 0.95, 0.85, 0.85, 0.95). However, it is **strongly recommended** to tune this parameter to any given dataset.}
  \item{use.sex.info}{Whether to process sex chromosome probes separately (recommended). If TRUE (default), uses sample sex information in determining how to call copy number with these probes. If the number of chrX/chrY probes per sample is low, male sex chromosome CNAs should be interpreted with caution. See \code{\link{process.xy.probes}} for more information.}
}
\details{Get copy number values for each probe using a reference of pooled normal samples. Preferable when most tumour samples do not have a matched normal. Pooled reference is created by taking the mean count of normal sample probes. Please note that results are re-centered around the neutral copy number (2 for autosomes and female sex chromosomes, and 1 for male sex chromosomes).}
\value{A list containing 4 items (each being a gene by sample matrix):
\cr1) copy number values rounded according to the method specified
\cr2) raw (not rounded) copy number values
\cr3) rounded copy number values for normal samples only
\cr4) raw copy number values for normal samples only}
\references{
  See NanoString website for PDFs on analysis guidelines:
    \url{http://www.nanostring.com/support/prod-lit/}

  The NanoString assay is described in the paper:
    Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\author{Emilie Lalonde and Dorota Sendorek}
\seealso{\code{\link{call.cnas.with.matched.normals}}}
