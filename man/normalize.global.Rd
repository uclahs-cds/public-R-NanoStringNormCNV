\name{normalize.global}
\alias{normalize.global}
\title{Wrapper function to call all normalization function requested}
\description{Wrapper function to call all normalization function requested using all data at once}
\usage{
normalize.global(raw.data, cc, bc, sc, do.nsn, do.rcc.inv, covs, plot.types = 'all')
}
\arguments{
  \item{raw.data}{The NanoString data-frame with raw counts for all samples in experiment}
  \item{cc}{Requested CodeCount NSN normalization. Can be 'none', 'sum', and 'geo.mean'.}
  \item{bc}{Requested Background NSN normalization. Can be 'none', 'mean', 'mean.2sd' and 'max'.}
  \item{sc}{Requested SampleContent NSN normalization. Can be 'none', 'housekeeping.geo.mean', 'total.sum', 'top.geo.mean' and 'low.cv.geo.mean'.}
  \item{do.nsn}{Whether or not to even run NanoStringNorm (NSN). Should be inferred from 3 previous parameters}
  \item{do.rcc.inv}{Whether or not to run invariant probe normalization. This is run separately from NSN, and is the recommended method by NanoString.}
  \item{covs}{A data-frame with covariates for which to assess batch effects in NanoStringNorm}
  \item{plot.types}{Which NSN plots to generate, if any, after normalization. Defaults to all}
}
\details{Normalize all data in study at once using NanoStringNorm and/or invariant probe normalizaiton. See individual normalization functions for more details}
\value{A data-frame with normalized data}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::NanoStringNorm}
