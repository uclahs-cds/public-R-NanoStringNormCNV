\name{make.cna.heatmap}
\alias{make.cna.heatmap}
\title{Make a heatmap for CNAs}
\description{Make a heatmap for CNAs}
\usage{
make.cna.heatmap(
	nano.cnas, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL,
  covs.legend = NULL,	rounded = FALSE, clust.dim = 'both', 
  centering.value = 2, ...
	)
}
\arguments{
  \item{nano.cnas}{A matrix with CNA info, genes X samples}
  \item{fname.stem}{To use in filename, description of how CNAs were called is suggested. Defaults to NULL}
  \item{covs.rows}{A covariate object for samples to be passed directly to create.heatmap. Defaults to NULL}
  \item{covs.cols}{A covariate object for genes to be passed directly to create.heatmap. Defaults to NULL}
  \item{covs.legend}{A legend object to be passed directly to create.heatmap. Defaults to NULL}
  \item{rounded}{Whether the CNA data has been rounded to the nearest integer. If false (default), correlation is used as the method. If true, jaccard is used}
  \item{clust.dim}{Which dimensions to cluster. Defaults to both}
  \item{centering.value}{The center value of the colour map. Defaults to 2}
  \item{...}{Other parameters to be passed to create.heatmap}
}
\details{Make a clustered heatmap for all samples and genes}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
