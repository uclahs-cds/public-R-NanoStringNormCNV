\name{invariant.probe.norm}
\alias{invariant.probe.norm}
\title{Normalize data using invariant probes}
\description{Normalize data using invariant probes to control for differences in sample input}
\usage{
invariant.probe.norm(nano.df)
}
\arguments{
  \item{nano.df}{The NanoString data-frame with normalized counts}
}
\details{Normalize data using invariant probes to control for differences in sample input. Counts are multiplied by a normalization factor (mean INV probe counts in sample / mean(mean INV probe count of each sample)).
}
\value{Nano.df with normalized counts}
\references{Manual pp 8-10}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::sample.content.normalization}
