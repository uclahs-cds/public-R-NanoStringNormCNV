\name{process.xy.probes}
\alias{process.xy.probes}
\title{Processes sex chromosome probes}
\description{Identifies and processes sex chromosome probe information}
\usage{
process.xy.probes(normalized.data, sex.info)
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts for all samples in experiment. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Sex chromosome probes are identified as probes whose names ('Name' column) contain either "chrX" or "chrY".}
  \item{sex.info}{A data-frame containing columns 'SampleID' and 'Sex' for each sample in 'normalized.data'. See output of \code{\link{load.phenodata}}}
}
\details{Sex chromosome probes are identified and processed according to user-provided sample sex information. Chromosome Y probes only are removed from female samples. Chromosome X and Y probes are removed from samples where sex information is not available. Finally, chromosome X and Y probe values for male samples are extracted and placed into a new data-frame, enabling separate downstream processing.}
\value{A list containing three items: 1) a data-frame holding 'normalized.data' with chrY probes removed from female samples, and chrX and chrY probes removed from both male samples and samples where sex is not available, 2) a data-frame holding 'normalized.data' for male chrX and chrY probes only (set to NULL if there are no male samples) and 3) names of all probes identified as corresponding to sex chromosomes}
\author{Dorota Sendorek}
