\name{call.cnas.with.pooled.normals}
\alias{call.cnas.with.pooled.normals}
\title{Call CNAs using a pooled normal reference}
\description{Perform CNA calling on both tumour and normal samples against a reference comprised of pooled normal samples.}
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
  \item{call.method}{The CNA calling method used to round tumour:normal ratios to copy number thresholds. Options are:
  \cr- 1 to use pre-defined thresholds derived from NanoString recommendations: 0.4, 1.5, 2.5, 3.5 (Default)
  \cr- 2 to use the kernel density method to calculate cohort-specific copy number thresholds for CNA calling (see 'kd.values' argument)
  \cr- 3 to use thresholds calculated from the reference sample set's median min and max values: min, min + SD, max - SD, max
  }
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values (between 0 and 1) to be used to calculate copy number thresholds (when 'call.method' 2 is requested). Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain. For probes identified as belonging to male sex chromosomes (neutral CN = 1), deletions will be determined using the heterozygous threshold. Default values are: 0.85, 0.95 (alternatively, this could be written as: 0.95, 0.85, 0.85, 0.95). However, it is **strongly recommended** to tune this parameter to any given dataset.}
  \item{use.sex.info}{Whether to process sex chromosome probes separately (recommended). If TRUE (default), uses sample sex information in determining how to call copy number with these probes. For female samples, chromosome Y probes are removed completely. For samples where sex information is not available, both chromosome X and Y probes are removed. If the number of chrX and chrY probes per sample is low, male sex chromosome CNAs should be interpreted with caution.}
}
\details{
  Get copy number values for each probe using a reference of pooled normal samples. Preferable when most tumour samples do not have a matched normal. Pooled reference is created by taking the mean count of normal sample probes. Please note that results are re-centered around the neutral copy number (2 for autosomes and female sex chromosomes, and 1 for male sex chromosomes).
}
\value{
  A list containing 4 items (each being a gene by sample matrix):
  \cr1) copy number values rounded according to the method specified
  \cr2) raw (not rounded) copy number values
  \cr3) rounded copy number values for normal samples only
  \cr4) raw copy number values for normal samples only
}
\references{
  See NanoString website for PDFs on analysis guidelines:
    \url{https://www.nanostring.com/support/product-support/support-documentation}

  The NanoString assay is described in the paper:
    Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\author{Emilie Lalonde and Dorota Sendorek}
\seealso{\code{\link{call.cnas.with.matched.normals}}}
\examples{\dontrun{
  # load data
  data(NanoString.DNA.norm);
  data(PhenoData);

  # call CNAs by rounding ratios to thresholds derived from NanoString recommendations
  cnas <- call.cnas.with.pooled.normals(
    normalized.data = NanoString.DNA.norm,
    phenodata = PhenoData,
    call.method = 1,
    use.sex.info = TRUE    
    );

  # call CNAs using kernel density method
  cnas <- call.cnas.with.pooled.normals(
    normalized.data = NanoString.DNA.norm,
    phenodata = PhenoData,
    call.method = 2,
    kd.values = c(0.99, 0.87, 0.89, 0.96),
    use.sex.info = TRUE
    );

  # call CNAs by rounding to thresholds calculated from the copy number range of 
  # the reference samples
  cnas <- call.cnas.with.pooled.normals(
    normalized.data = NanoString.DNA.norm,
    phenodata = PhenoData,
    call.method = 3,
    use.sex.info = TRUE
    );  
}}
