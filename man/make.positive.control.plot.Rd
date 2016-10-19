\name{make.positive.control.plot}
\alias{make.positive.control.plot}
\title{Make heatmaps of correlations}
\description{Make heatmaps of correlations between target concentration and actual counts}
\usage{
make.positive.control.plot(correlations, covs = NULL, print.x.labels = TRUE)
}
\arguments{
  \item{correlations}{A data-frame with columns for sample IDs ('SampleID') and Pearson correlations between target concentration and actual counts ('R2'). See output of \code{\link{positive.control.norm}}}
  \item{covs}{A sample by covariate data-frame. Currently only accepts sample covariates 'type' and 'cartridge'. A 'SampleID' column is mandatory to match to 'correlations' sample IDs. Defaults to NULL}
  \item{print.x.labels}{Whether to display the sample IDs on the heatmap's x-axis. Defaults to TRUE}
}
\details{Make heatmaps of correlations between target concentration and actual counts. First heatmap shows full range of correlation, and second is zoomed in to observed values only}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
