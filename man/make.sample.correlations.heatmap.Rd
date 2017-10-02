\name{make.sample.correlations.heatmap}
\alias{make.sample.correlations.heatmap}
\title{Make a heatmap of the inter-sample correlations}
\description{Make a heatmap of inter-sample correlations of code counts}
\usage{
make.sample.correlations.heatmap(
  nano.counts, cor.method = 'pearson', fname.stem = NULL, covs = NULL
  )
}
\arguments{
  \item{nano.counts}{A gene by sample matrix or data-frame of NanoString counts, where row names are probe names and column names are sample IDs. Counts can be normalized or raw.}
  \item{cor.method}{Which distance method to use for clustering. Accepts 'pearson' (default), 'kendall', or 'spearman'}
  \item{fname.stem}{To use in filename. Description of how CNAs were called is suggested. Defaults to NULL}
  \item{covs}{A sample by covariate data-frame. Currently only accepts sample covariates 'Type' and 'Cartridge' (see output of \code{\link{load.phenodata}}). A 'SampleID' column is mandatory to match to 'nano.counts' sample IDs. Defaults to NULL}
}
\details{Make a clustered heatmap of code count correlations for all samples}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
\examples{\dontrun{
  # load data
  data(NanoString.DNA.raw);
  data(PhenoData);

	# make sure raw count data frame has gene names for rownames!
	NanoString.DNA.formatted <- NanoString.DNA.raw[, -(1:3)];
	rownames(NanoString.DNA.formatted) <- NanoString.DNA.raw$Name;

	# plot raw NanoString count correlations
	make.sample.correlations.heatmap(
		nano.counts = NanoString.DNA.formatted,
		covs = PhenoData[, c('SampleID', 'Cartridge', 'Type')]
		);
}}
