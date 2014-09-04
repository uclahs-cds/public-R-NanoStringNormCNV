\name{make.counts.heatmap}
\alias{make.counts.heatmap}
\title{Make a heatmap for NanoString counts}
\description{Make a heatmap for counts}
\usage{
make.counts.heatmap(nano.counts, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL)
}
\arguments{
  \item{nano.counts}{A matrix with nanostring count info, genes X samples}
  \item{fname.stem}{To use in filename, description of how CNAs were called is suggested}
  \item{covs.rows}{A covariate object for samples to be passed directly to create.heatmap}
  \item{covs.cols}{A covariate object for genes to be passed directly to create.heatmap}
  \item{covs.legend}{A legend object to be passed directly to create.heatmap}
}
\details{Make a clustered heatmap for all samples and genes}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
