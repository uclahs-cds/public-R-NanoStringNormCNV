\name{normalize.global}
\alias{normalize.global}
\title{Wrapper function to normalize raw data, globally}
\description{Wrapper function to apply all normalization functions requested to samples globally}
\usage{
  normalize.global(
    raw.data, cc, bc, sc, oth, do.rcc.inv, covs,
    transform.data = TRUE, plot.types = 'all', phenodata = NULL
    )
}
\arguments{
  \item{raw.data}{A gene by sample data-frame of raw NanoString counts for all samples in experiment. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{cc}{Requested CodeCount NSN normalization. Can be 'none', 'sum', and 'geo.mean'.}
  \item{bc}{Requested Background NSN normalization. Can be 'none', 'mean', 'mean.2sd' and 'max'.}
  \item{sc}{Requested SampleContent NSN normalization. Can be 'none', 'housekeeping.geo.mean', 'total.sum', 'top.geo.mean' and 'low.cv.geo.mean'. To run 'housekeeping.geo.mean', user must ensure the NanoString dataset includes probes belonging to the 'Housekeeping' code class. These can be probes that were found to have consistently low variance across data.}
  \item{oth}{Requested OtherNorm NSN normalization. Can be 'none', 'vsn', 'quant' and 'rank'.}
  \item{do.rcc.inv}{Whether or not to run invariant probe normalization. This method corrects for sample content and can be used in place of 'sc'. It is run separately from NSN and is the technique recommended in NanoString guidelines.}
  \item{covs}{A data-frame of covariates with which to assess batch effects in NanoStringNorm. Passed directly to NanoStringNorm::NanoStringNorm 'traits' argument (see function for specifics). Set to NA if none.}
  \item{transform.data}{Whether to transform data so output does not contain negative values. Defaults to TRUE}
  \item{plot.types}{Which NSN plots to generate, if any, after normalization. Passed directly to NanoStringNorm::Plot.NanoStringNorm. Defaults to 'all'}
  \item{phenodata}{A data-frame of sample annotation which includes the fields 'SampleID' and 'Type'. Recommended to use output of \code{\link{load.phenodata}}. Defaults to NULL}
}
\details{
  Normalizes all samples in experiment together using NanoStringNorm and/or invariant probe normalizaiton. Normalization methods requested are implemented in the following order: code count, background correction, sample content, invariant probe, and other. See NanoStringNorm and \code{\link{invariant.probe.norm}} for more details.
}
\value{A gene by sample data-frame of normalized counts, where first three columns contain 'CodeClass', 'Name', and 'Accession' information}
\references{
  See NanoString website for PDFs on analysis guidelines:
    \url{https://www.nanostring.com/support/product-support/support-documentation}

  The NanoString assay is described in the paper:
    Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).

  The NanoStringNorm package is described in the paper:
    Waggott, D., Chu, K., Yin, S., Wouters, B.G., Liu, F.F. & Boutros, P.C. NanoStringNorm: an extensible R package for the pre-processing of NanoString mRNA and miRNA data. Bioinformatics 28(11), 1546-1548 (2012).
}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::NanoStringNorm, \code{\link{normalize.per.chip}}}
\examples{\dontrun{
  # load data
  data(NanoString);
  data(PhenoData);

  # perform only invariant probe normalization
  NanoString.norm <- normalize.global(
    raw.data = NanoString,
    cc = 'none',
    bc = 'none',
    sc = 'none',
    oth = 'none',
    do.rcc.inv = TRUE,
    covs = NA,
    phenodata = PhenoData
    );

  ## accounting for batch effects
  # include covariates for sample cartridge and sample type 
  # must be binary as these are passed directly to NanoStringNorm 'traits'
  covs <- as.data.frame(matrix(
    1,
    nrow = nrow(PhenoData),
    ncol = length(unique(PhenoData$Cartridge)),
    dimnames = list(
      PhenoData$SampleID,
      paste0("Cartridge", unique(PhenoData$Cartridge))
      )
    ));

  for (n in 1:nrow(PhenoData)) {
    covs[n, which(unique(PhenoData$Cartridge) == PhenoData$Cartridge[n])] <- 2;
    }

  covs$Type <- ifelse(PhenoData$Type == 'Reference', 1, 2);

  # normalize
  NanoString.norm <- normalize.global(
    raw.data = NanoString,
    cc = 'none',
    bc = 'none',
    sc = 'none',
    oth = 'none',
    do.rcc.inv = TRUE,
    covs = covs,
    phenodata = PhenoData
    );
}}
