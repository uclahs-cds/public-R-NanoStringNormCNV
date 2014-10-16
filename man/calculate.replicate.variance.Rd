\name{calculate.replicate.variance}
\alias{calculate.replicate.variance}
\title{Determine the variance of replicate samples}
\description{Return a variance matrix for replicate samples for all genes}
\usage{
calculate.replicate.variance(input, pheno, gene.names, var.function = 'sd')
}
\arguments{
  \item{input}{The NanoString normalized counts for all replicate samples}
  \item{pheno}{A data-frame containing the Name and SampleID for each replicate sample. The column names of input should be the SampleIDs (there will be duplicate column names)}
  \item{gene.names}{A vector with the names of the genes ordered as in the input variable.}
  \item{var.function}{Which variance-type function should be used? Defaults to sd}
}
\details{Determine the variance for each gene  within replicate samples.}
\value{A matrix with n_replicates columns and n_genes rows, and each cell is the variance of that gene in that sample's replicates}
\author{Cindy Yao and Emilie Lalonde}
\seealso{calculate.replicate.variance}
