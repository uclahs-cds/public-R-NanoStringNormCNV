\name{make.sample.correlations.heatmap}
\alias{make.sample.correlations.heatmap}
\title{Make a heatmap of the inter-sample correlations}
\description{Make a heatmap for inter-sample correlations based on code counts}
\usage{
make.sample.correlations.heatmap(
  nano.counts, cor.method = 'pearson', fname.stem = NULL, covs = NULL
  )
}
\arguments{
  \item{nano.counts}{A gene by sample matrix or data-frame of NanoString counts, where row names are probe names and column names are sample IDs. Counts can be normalized or raw.}
  \item{cor.method}{Which distance method to use for clustering. Accepts 'pearson' (default), 'kendall', or 'spearman'}
  \item{fname.stem}{To use in filename, description of data is suggested. Defaults to NULL}
  \item{covs}{A sample by covariate data-frame. Currently only accepts sample covariates 'Type' and 'Cartridge'. A 'SampleID' column is mandatory to match to 'nano.counts' sample IDs. Defaults to NULL}
}
\details{Make a clustered heatmap for all samples based on code count data}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
