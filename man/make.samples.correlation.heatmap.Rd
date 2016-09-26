\name{make.samples.correlation.heatmap}
\alias{make.samples.correlation.heatmap}
\title{Make a heatmap of the inter-sample correlations}
\description{Make a heatmap for inter sample correlations based on code counts}
\usage{
make.samples.correlation.heatmap(data, cor.method = 'pearson', fname, covs, covs.legend)
}
\arguments{
  \item{data}{A matrix with code count info, genes X samples. Can be normalized or raw.}
  \item{cor.method}{Which distance method to use for clustering. Defaults to 'pearson'}
  \item{fname}{The filename to use}
  \item{covs}{A covariate object for samples to be passed directly to create.heatmap}
  \item{covs.legend}{A legend object to be passed directly to create.heatmap}
}
\details{Make a clustered heatmap for all samples  based on code count data}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
