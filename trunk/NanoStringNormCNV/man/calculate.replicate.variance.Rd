\name{calculate.replicate.variance}
\alias{calculate.replicate.variance}
\title{Determine the variance of replicate samples}
\description{Return a variance matrix for replicate samples for all genes}
\usage{
calculate.replicate.variance(norm.data.reps, phenodata.reps, var.function = 'sd')
}
\arguments{
  \item{norm.data.reps}{A gene by sample data-frame of normalized NanoString counts for each replicate sample. Column names must match sample IDs found in 'phenodata.reps'. Expected row names are gene probe names}
  \item{phenodata.reps}{A data-frame containing the 'Name' and 'SampleID' for each replicate sample}
  \item{var.function}{Which variance-type function should be used. Defaults to sd}
}
\details{Determine the variance for each gene within replicate samples.}
\value{A matrix with n_replicates columns and n_genes rows where each cell is the variance of given gene in given sample's replicates}
\author{Cindy Yao and Emilie Lalonde}
\seealso{calculate.replicate.variance}
