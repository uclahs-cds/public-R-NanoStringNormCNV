\name{make.cna.heatmap}
\alias{make.cna.heatmap}
\title{Make a heatmap of CNAs}
\description{Make a heatmap of CNAs}
\usage{
make.cna.heatmap(
  nano.cnas, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL,
  rounded = FALSE, clust.dim = 'both', centering.value = 2, min.cn = NULL, ...
  )
}
\arguments{
  \item{nano.cnas}{A gene by sample matrix of CNA information (raw ratios or rounded CNA calls)}
  \item{fname.stem}{To use in filename. Description of how CNAs were called is suggested. Defaults to NULL}
  \item{covs.rows}{A sample by covariate data-frame. Currently only accepts sample covariates 'Type' and 'Cartridge'. A 'SampleID' column is mandatory to match to 'nano.counts' sample IDs. Defaults to NULL}
  \item{covs.cols}{A gene by covariate data-frame. Currently only accepts gene covariate 'CodeClass'. A 'Name' column is mandatory to match to 'nano.counts' gene probe names. Defaults to NULL}
  \item{rounded}{Whether the CNA data has been rounded to the nearest integer. If false (default), euclidean distance measure is used to cluster data. If true, jaccard is used.}
  \item{clust.dim}{Which dimensions to cluster. Accepts 'none', 'rows', 'columns', or 'both' (default)}
  \item{centering.value}{The center value of the colour map. Defaults to 2}
  \item{min.cn}{A value to be set as the colour scheme minimum. If NULL (default), takes minimum value of 'nano.cnas'}
  \item{...}{Other parameters to be passed to the BoutrosLab.plotting.general::create.heatmap function}
}
\details{Make a clustered heatmap of CNAs for all samples and genes}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
