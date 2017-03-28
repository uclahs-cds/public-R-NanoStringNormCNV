\name{make.cna.densities.plots}
\alias{make.cna.densities.plots}
\title{Make density plots for normalized counts}
\description{Make density plots of copy number calls per sample and per gene}
\usage{
make.cna.densities.plots(nano.cnas, fname.stem = NULL)
}
\arguments{
  \item{nano.cnas}{A gene by sample matrix of CNA information (raw ratios or rounded CNA calls)}
  \item{fname.stem}{To use in filename. Description of how CNAs were called is suggested. Defaults to NULL}
}
\details{Make density plots of copy number calls. Done separately for all genes and separately for all samples}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
