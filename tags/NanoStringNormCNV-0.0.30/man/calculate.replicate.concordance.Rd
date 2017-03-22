\name{calculate.replicate.concordance}
\alias{calculate.replicate.concordance}
\title{Determine concordance across replicate samples}
\description{Quantifies CNA call concordance for all genes/probes found across replicate samples}
\usage{
calculate.replicate.concordance(cna.rounded.reps, phenodata.reps)
}
\arguments{
  \item{cna.rounded.reps}{A gene by sample data-frame of CNA calls for all replicate samples. Column names must match sample IDs found in 'phenodata.reps'.}
  \item{phenodata.reps}{A data-frame containing phenotype information for all replicate samples. Must include columns 'Name' and 'SampleID'. Replicates are identified by their matching sample names (phenodata column 'Name'). See output of \code{\link{load.phenodata}}}
}
\value{A gene by sample name data-frame of concordance values (1 signifies matching CN values whereas 0 signifies discordant CN values)}
\author{Cindy Yao and Emilie Lalonde}
\seealso{\code{\link{calculate.replicate.variance}}}
