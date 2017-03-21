\name{calculate.replicate.variance}
\alias{calculate.replicate.variance}
\title{Determine variance across replicate samples}
\description{Return a normalized counts variance matrix for all genes/probes found in replicate samples}
\usage{
calculate.replicate.variance(normalized.data.reps, phenodata.reps, var.function = 'sd')
}
\arguments{
  \item{normalized.data.reps}{A gene by sample data-frame of normalized NanoString counts for all replicate samples. Column names must match sample IDs found in 'phenodata.reps'.}
  \item{phenodata.reps}{A data-frame containing phenotype information for all replicate samples. Must include columns 'Name' and 'SampleID'. See output of \code{\link{load.phenodata}}}
  \item{var.function}{Which variance-type function should be used. Defaults to 'sd'}
}
\details{Determine the variance of normalized counts for each gene across replicate samples.}
\value{A gene by sample name data-frame of variances}
\author{Cindy Yao and Emilie Lalonde}
\seealso{\code{\link{calculate.replicate.variance}}}
