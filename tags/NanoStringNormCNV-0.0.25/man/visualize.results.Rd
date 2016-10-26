\name{visualize.results}
\alias{visualize.results}
\title{visualize results}
\description{Wrapper function that creates a variety of plots to display results}
\usage{
visualize.results(
  raw.counts, norm.counts, phenodata = NULL, cna.rounded = NULL,
  cna.raw = NULL, max.cn = 5, replicate.eval = NULL, exclude.covs = FALSE
  )
}
\arguments{
  \item{raw.counts}{A gene by sample data-frame of raw NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{norm.counts}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{Data-frame containing annotation for sample covariates. Expects columns 'SampleID' and 'type' or 'cartridge' or both. Defaults to NULL}
  \item{cna.rounded}{A gene by sample data-frame of rounded CNA calls, where column names are sample IDs and row names are gene probe names. Defaults to NULL}
  \item{cna.raw}{A gene by sample data-frame of raw CNA calls, where column names are sample IDs and row names are gene probe names. Defaults to NULL}
  \item{max.cn}{Maximum copy number for 'cna.raw' data. Sets copy numbers of more than 'max.cn' to 'max.cn'. Defaults to 5}
  \item{replicate.eval}{A list holding results of replicate evaluation. See output of \code{\link{evaluate.replicates}}}
  \item{exclude.covs}{Whether to exclude covariates from plots. Defaults to FALSE}
}
\details{Plots are saved to files in current working directory.}
\value{None}
\author{Emilie Lalonde and Dorota Sendorek}
