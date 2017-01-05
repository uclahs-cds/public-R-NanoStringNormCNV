\name{calculate.replicate.concordance}
\alias{calculate.replicate.concordance}
\title{Determine concordance across replicate samples}
\description{Return a CNA call concordance matrix for all genes/probes found in replicate samples}
\usage{
calculate.replicate.concordance(cna.rounded.reps, phenodata.reps)
}
\arguments{
  \item{cna.rounded.reps}{A gene by sample data-frame of CNA calls for all replicate samples. Column names must match sample IDs found in 'phenodata.reps'.}
  \item{phenodata.reps}{A data-frame containing phenotype information for all replicate samples. Must include columns 'Name' and 'SampleID'. See output of \code{\link{load.phenodata}}}
}
\details{Determine whether CNA calls for each gene are concordant between replicate samples.}
\value{A gene by sample name data-frame where each cell contains either 1 for concordance or 0 for discordance}
\author{Cindy Yao and Emilie Lalonde}
\seealso{calculate.replicate.variance}
