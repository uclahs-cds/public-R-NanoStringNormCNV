\name{calculate.replicate.concordance}
\alias{calculate.replicate.concordance}
\title{Determine the concordance of replicate samples}
\description{Return a concordance matrix for replicate samples for all genes}
\usage{
calculate.replicate.concordance(input, pheno, gene.names)
}
\arguments{
  \item{input}{The NanoString CNA calls for all replicate samples}
  \item{pheno}{A data-frame containing the Name and SampleID for each replicate sample. The column names of input should be the SampleIDs (there will be duplicate column names)}
  \item{gene.names}{A vector with the names of the genes ordered as in the input variable.}
}
\details{Determine whether each gene is concordant within replicate samples.}
\value{A matrix with n_replicates columns and n_genes rows, and each cell is 1 for concordance or 0 for discordance}
\author{Cindy Yao and Emilie Lalonde}
\seealso{calculate.replicate.variance}
