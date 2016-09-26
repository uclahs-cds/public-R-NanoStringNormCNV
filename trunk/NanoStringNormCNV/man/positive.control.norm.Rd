\name{positive.control.norm}
\alias{positive.control.norm}
\title{Do QAQC based on positive controls}
\details{Find the correlation of positive control probes and the expected amounts}
\usage{
positive.control.norm(nano.df)
}
\arguments{
  \item{nano.df}{A data-frame with the raw counts}
}
\description{Find the correlation of positive control probes and the expected amounts based on NanoString manual (y=171.23x + 214.12)}
\value{A data-frame with SampleName and R2 for Pearson correlations}
\author{Cindy Yao and Emilie Lalonde}
