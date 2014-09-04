\name{make.positive.control.plot}
\alias{make.positive.control.plot}
\title{Make heatmaps of correlations}
\description{Make heatmaps of correlation between target concentration and actual counts}
\usage{
make.positive.control.plot(corr.df, fname.stem)
}
\arguments{
  \item{corr.df}{The correlation data-frame with the correlation in the first column}
  \item{fname.stem}{Base name for files}
}
\details{Make heatmaps of correlation between target concentration and actual counts. First heatmap shows full range of correlation, and second is zoomed in to observed values only}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
