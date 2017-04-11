\name{visualize.results}
\alias{visualize.results}
\title{visualize results}
\description{Wrapper function that creates a variety of plots, displaying results}
\usage{
visualize.results(
  raw.data, normalized.data, phenodata = NULL, cna.rounded = NULL,
  cna.raw = NULL, max.cn = 5, replicate.eval = NULL, exclude.covs = FALSE
  )
}
\arguments{
  \item{raw.data}{A gene by sample data-frame of raw NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{A data-frame of sample annotation. Must include columns 'SampleID' and 'Type' and/or 'Cartridge'. Recommended to use output of \code{\link{load.phenodata}}. Defaults to NULL}
  \item{cna.rounded}{A gene by sample data-frame of CNA calls. See 'cna.rounded' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL}
  \item{cna.raw}{A gene by sample data-frame of tumour-normal count ratios. See 'cna.raw' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL}
  \item{max.cn}{Maximum allowable ratio in 'cna.raw' data. Sets copy numbers of more than 'max.cn' to 'max.cn'. Defaults to 5}
  \item{replicate.eval}{A list of replicate evaluation results. Use output of \code{\link{evaluate.replicates}}}
  \item{exclude.covs}{Whether to exclude covariate information (provided in 'phenodata') from plots. Defaults to FALSE}
}
\details{Plots are saved in tiff files to working directory}
\value{None}
\author{Emilie Lalonde and Dorota Sendorek}
