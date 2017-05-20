\name{call.cnas.with.matched.normals}
\alias{call.cnas.with.matched.normals}
\title{Call CNAs using matched normal references}
\description{Perform CNA calling on tumour and matched normal samples only}
\usage{
call.cnas.with.matched.normals(
  normalized.data, phenodata, per.chip = FALSE, call.method = 1,
  kd.values = c(0.85, 0.95), use.sex.info = TRUE
  )
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probes are identified by "chrX" and "chrY" labels in their 'Name' and/or 'Accession' values.}
  \item{phenodata}{A data-frame of sample annotation. Recommended to use output of \code{\link{load.phenodata}}}
  \item{per.chip}{Whether to use reference samples specific to each chip when calculating the tumour:normal ratios from which copy numbers are generated. If TRUE, will skip CNA calling for matched tumour-normal samples that were processed on different chips. Defaults to FALSE.}
  \item{call.method}{The CNA calling method used to round tumour:normal ratios to copy number thresholds. Options are:
  \cr- 1 to use pre-defined thresholds derived from NanoString recommendations: 0.4, 1.5, 2.5, 3.5 (Default)
  \cr- 2 to use the kernel density method to calculate sample-specific copy number thresholds for CNA calling (see 'kd.values' argument)
  }
  \item{kd.values}{A 2- or 4-element vector containing the kernel density values (between 0 and 1) to be used to calculate copy number thresholds (when 'call.method' 2 is requested). Order is important. For 4 elements, values define homozygous deletion, heterozygous deletion, single copy gain, and multi-copy gain. For 2 elements, values define heterozygous deletion/gain and homozygous deletion/gain. For probes identified as belonging to male sex chromosomes (neutral CN = 1), deletions will be determined using the heterozygous threshold. Default values are: 0.85, 0.95 (alternatively, this could be written as: 0.95, 0.85, 0.85, 0.95). However, it is **strongly recommended** to tune this parameter to any given dataset.}
  \item{use.sex.info}{Whether to process sex chromosome probes separately (recommended). If TRUE (default), uses sample sex information in determining how to call copy number with these probes. For female samples, chromosome Y probes are removed completely. For samples where sex information is not available, both chromosome X and Y probes are removed. If the number of chrX and chrY probes per sample is low, male sex chromosome CNAs should be interpreted with caution.}
}
\details{Get copy number values for each probe using corresponding matched normal samples. Be aware that tumour samples without a matched normal get dropped.}
\value{A list containing 2 items: 1) a gene by sample data-frame of observed copy numbers rounded according to the method specified and 2) a gene by sample data-frame of raw (not rounded) observed copy numbers}
\references{
  See NanoString website for PDFs on analysis guidelines:
    \url{https://www.nanostring.com/support/product-support/support-documentation}

  The NanoString assay is described in the paper:
    Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}

\author{Emilie Lalonde and Dorota Sendorek}
\seealso{\code{\link{call.cnas.with.pooled.normals}}}
\examples{\dontrun{
  # load data
  data(NanoString.norm);
  data(PhenoData);

  # call CNAs by rounding ratios to thresholds derived from NanoString recommendations
  cnas <- call.cnas.with.matched.normals(
    normalized.data = NanoString.norm,
    phenodata = PhenoData,
    call.method = 1,
    use.sex.info = TRUE    
    );

  # call CNAs using kernel density method
  cnas <- call.cnas.with.matched.normals(
    normalized.data = NanoString.norm,
    phenodata = PhenoData,
    call.method = 2,
    kd.values = c(0.99, 0.87, 0.89, 0.96),
    use.sex.info = TRUE
    );
}}
