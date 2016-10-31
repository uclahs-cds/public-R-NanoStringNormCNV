\name{normalize.per.chip}
\alias{normalize.per.chip}
\title{Wrapper function to call all normalization function requested}
\description{Wrapper function to call all normalization function requested, processing one cartridge at a time}
\usage{
normalize.per.chip(
	pheno, raw.data, cc, bc, sc, oth, do.nsn, do.rcc.inv, covs,
  transform.data = TRUE, plot.types = 'all'
	)
}
\arguments{
  \item{pheno}{A data-frame with the annotations for each sample, which must include a 'cartridge' field. May include a 'type' field (specifying whether Tumour or Reference).}
  \item{raw.data}{A gene by sample data-frame of raw NanoString counts for all samples in experiment. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{cc}{Requested CodeCount NSN normalization. Can be 'none', 'sum', and 'geo.mean'.}
  \item{bc}{Requested Background NSN normalization. Can be 'none', 'mean', 'mean.2sd' and 'max'.}
  \item{sc}{Requested SampleContent NSN normalization. Can be 'none', 'housekeeping.geo.mean', 'total.sum', 'top.geo.mean' and 'low.cv.geo.mean'.}
  \item{oth}{Requested OtherNorm NSN normalization. Can be 'vsn', 'quant' and 'rank'.}
  \item{do.nsn}{Whether or not to even run NanoStringNorm (NSN). Should be inferred from 3 previous parameters}
  \item{do.rcc.inv}{Whether or not to run invariant probe normalization. This is run separately from NSN, and is the recommended method by NanoString.}
  \item{covs}{A data-frame with covariates for which to assess batch effects in NanoStringNorm}
  \item{transform.data}{Whether to transform data so output does not contain negative values. Defaults to TRUE}
  \item{plot.types}{Which NSN plots to generate, if any, after normalization. Passed directly to NanoStringNorm::Plot.NanoStringNorm. Defaults to 'all'}
}
\details{Normalize all data in study one cartridge at a time using NanoStringNorm and/or invariant probe normalizaiton. See individual normalization functions for more details.}
\value{A data-frame with normalized data, with all cartridges combined}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::NanoStringNorm}
