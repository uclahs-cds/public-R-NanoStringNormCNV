\name{calculate.replicate.variance}
\alias{calculate.replicate.variance}
\title{Determine count variance across replicate samples}
\description{Quantifies normalized count variance for all genes/probes found across sample replicates}
\usage{
calculate.replicate.variance(normalized.data.reps, phenodata.reps, var.function = 'sd')
}
\arguments{
  \item{normalized.data.reps}{A gene by sample data-frame of normalized NanoString counts for all replicate samples. Column names must match sample IDs found in 'phenodata.reps'.}
  \item{phenodata.reps}{A data-frame containing sample annotation for all replicate samples. Must include columns 'Name' and 'SampleID'. Replicates belonging to a single sample are identified as such by their matching sample names (phenodata column 'Name'). See output of \code{\link{load.phenodata}}}
  \item{var.function}{The variance-type function to be used. Defaults to 'sd'}
}
\value{A gene by sample name data-frame of variances}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNormCNV:::calculate.replicate.concordance}
\keyword{internal}
