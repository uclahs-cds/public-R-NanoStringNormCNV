\name{normalize.per.chip}
\alias{normalize.per.chip}
\title{Wrapper function to call all normalization function requested}
\description{Wrapper function to apply all normalization functions requested to data by cartridge}
\usage{
normalize.per.chip(
	phenodata, raw.data, cc, bc, sc, oth, do.nsn, do.rcc.inv, covs,
  transform.data = TRUE, plot.types = 'all'
	)
}
\arguments{
  \item{phenodata}{A data-frame with the annotations for each sample, which must include a 'Cartridge' field. May include a 'Type' field (specifying whether Tumour or Reference) for plotting invariant probe data.}
  \item{raw.data}{A gene by sample data-frame of raw NanoString counts for all samples in experiment. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{cc}{Requested CodeCount NSN normalization. Can be 'none', 'sum', and 'geo.mean'.}
  \item{bc}{Requested Background NSN normalization. Can be 'none', 'mean', 'mean.2sd' and 'max'.}
  \item{sc}{Requested SampleContent NSN normalization. Can be 'none', 'housekeeping.geo.mean', 'total.sum', 'top.geo.mean' and 'low.cv.geo.mean'.}
  \item{oth}{Requested OtherNorm NSN normalization. Can be 'vsn', 'quant' and 'rank'.}
  \item{do.nsn}{Whether or not to even run NanoStringNorm (NSN). Should be inferred from 3 previous parameters}
  \item{do.rcc.inv}{Whether or not to run invariant probe normalization. This is run separately from NSN, and is the recommended method by NanoString.}
  \item{covs}{A data-frame with covariates for which to assess batch effects in NanoStringNorm. Passed directly to NanoStringNorm::NanoStringNorm 'traits' argument, see function for more information.}
  \item{transform.data}{Whether to transform data so output does not contain negative values. Defaults to TRUE}
  \item{plot.types}{Which NSN plots to generate, if any, after normalization. Passed directly to NanoStringNorm::Plot.NanoStringNorm. Defaults to 'all'}
}
\details{
  Normalizes all data in study one cartridge at a time using NanoStringNorm and/or invariant probe normalization. Normalization methods requested are performed in the following order: code count, background correction, sample content, invariant probe, and other. See individual normalization functions for details.
}
\value{A gene by sample data-frame of normalized counts, where first three columns contain 'CodeClass', 'Name', and 'Accession' information}
\references{
  See NanoString website for PDFs on analysis guidelines:
    \url{http://www.nanostring.com/support/prod-lit/}

  The NanoString assay is described in the paper:
    Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).

  The NanoStringNorm package is described in the paper:
    Waggott, D., Chu, K., Yin, S., Wouters, B.G., Liu, F.F. & Boutros, P.C. NanoStringNorm: an extensible R package for the pre-processing of NanoString mRNA and miRNA data. Bioinformatics 28(11), 1546-1548 (2012).
}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::NanoStringNorm}
