\name{process.xy.probes}
\alias{process.xy.probes}
\title{Processes sex chromosome probes}
\description{Identifies and processes sex chromosome probe information prior to CNA calling}
\usage{
process.xy.probes(normalized.data, sex.info)
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts for all samples in experiment. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probes are identified as those containing "chrX" and "chrY" labels in their 'Name' and/or 'Accession' values.}
  \item{sex.info}{A data-frame containing columns 'SampleID' and 'Sex' for each sample in 'normalized.data'. Recommended to use output of \code{\link{load.phenodata}}}
}
\details{Sex chromosome probes are processed according to user-provided sample sex information. For female samples, chromosome Y probes are removed completely. For samples where sex information is not available, both chromosome X and Y probes are removed. For male samples, chromosome X and Y probe values are extracted and placed into a new data-frame, enabling separate downstream analysis.}
\value{A list containing three items: 1) a data-frame holding 'normalized.data' with chrY probes removed from female samples, and chrX and chrY probes removed from both male samples and samples where sex is not available, 2) a data-frame holding 'normalized.data' for male chrX and chrY probes only (set to NULL if there are no male samples) and 3) names of all probes identified as belonging to sex chromosomes}
\author{Dorota Sendorek}
