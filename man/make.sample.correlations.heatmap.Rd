\name{make.sample.correlations.heatmap}
\alias{make.sample.correlations.heatmap}
\title{Make a heatmap of the inter-sample correlations}
\description{Make a heatmap for inter sample correlations based on code counts}
\usage{
make.sample.correlations.heatmap(
	data, cor.method = 'pearson', fname.stem = NULL,
	covs = NULL, covs.legend = NULL
	)
}
\arguments{
  \item{data}{A matrix with code count info, genes X samples. Can be normalized or raw.}
  \item{cor.method}{Which distance method to use for clustering. Defaults to 'pearson'}
  \item{fname.stem}{To use in filename, description of data is suggested. Defaults to NULL}
  \item{covs}{A covariate object for samples to be passed directly to create.heatmap. Defaults to NULL}
  \item{covs.legend}{A legend object to be passed directly to create.heatmap. Defaults to NULL}
}
\details{Make a clustered heatmap for all samples  based on code count data}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
