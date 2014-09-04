\name{make.cna.heatmap}
\alias{make.cna.heatmap}
\title{Make a heatmap for CNAs}
\description{Make a heatmap for CNAs}
\usage{
make.cna.heatmap(nano.cnas, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL)
}
\arguments{
  \item{nano.cnas}{A matrix with CNA info, genes X samples}
  \item{fname.stem}{To use in filename, description of how CNAs were called is suggested}
  \item{covs.rows}{A covariate object for samples to be passed directly to create.heatmap}
  \item{covs.cols}{A covariate object for genes to be passed directly to create.heatmap}
  \item{covs.legend}{A legend object to be passed directly to create.heatmap}
}
\details{Make a clustered heatmap for all samples and genes}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
