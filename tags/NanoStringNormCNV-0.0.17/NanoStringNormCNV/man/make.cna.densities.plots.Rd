\name{make.cna.densities.plots}
\alias{make.cna.densities.plots}
\title{Make density plots for normalized counts}
\description{Make density plots for normalized CNA counts per sample and per gene}
\usage{
make.cna.densities.plots(cnas, fname.stem, xlab)
}
\arguments{
  \item{cnas}{A matrix with CNA info, genes X samples}
  \item{fname.stem}{To use in filename, description of how CNAs were called is suggested}
  \item{xlab}{A description of type of CNA data}
}
\details{Make density plots of normalized counts after reference subtraction, once for all genes and once for all samples}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
