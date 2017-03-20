\name{make.positive.control.plot}
\alias{make.positive.control.plot}
\title{Make heatmaps of correlations}
\description{Make heatmaps of correlations between target concentrations and actual counts of positive control probes.}
\usage{
make.positive.control.plot(correlations, covs = NULL, print.x.labels = TRUE)
}
\arguments{
  \item{correlations}{A data-fram of Pearson R squared values by sample ID. See output of \code{\link{positive.control.qc}}}
  \item{covs}{A sample by covariate data-frame. Currently only accepts sample covariates 'Type' and 'Cartridge'. A 'SampleID' column is mandatory to match to sample IDs in 'correlations'. Defaults to NULL}
  \item{print.x.labels}{Whether to display the sample IDs on the heatmap's x-axis. Defaults to TRUE}
}
\details{Creates two heatmaps that display the sample R squared values, ascending.}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
