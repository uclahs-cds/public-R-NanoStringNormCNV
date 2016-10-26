\name{calculate.replicate.concordance}
\alias{calculate.replicate.concordance}
\title{Determine the concordance of replicate samples}
\description{Return a concordance matrix for replicate samples for all genes}
\usage{
calculate.replicate.concordance(cnas.reps, phenodata.reps)
}
\arguments{
  \item{cnas.reps}{A gene by sample data-frame of NanoString CNA calls for each replicate sample. Column names must match sample IDs found in 'phenodata.reps'. Expected row names are gene probe names}
  \item{phenodata.reps}{A data-frame containing the Name and SampleID for each replicate sample}
}
\details{Determine whether each gene is concordant within replicate samples.}
\value{A matrix with n_replicates columns and n_genes rows where each cell contains either 1 for concordance or 0 for discordance}
\author{Cindy Yao and Emilie Lalonde}
\seealso{calculate.replicate.variance}
