\name{get.sample.specific.cna.thresholds}
\alias{get.sample.specific.cna.thresholds}
\title{Determine log2-ratio thresholds to call CNAs}
\description{Determines sample-specific thresholds to call CNA states. Basic formula is median - x*SD.}
\usage{get.sample.specific.cna.thresholds(cna.data, percent = 0.8)}
\arguments{
  \item{cna.data}{A vector with the log2ratio values, no NAs allowed}
  \item{percent}{The percent of density you would like to have between the threshold bounds (0.1 - 0.99), Default = 0.8}
}
\details{This script will return the log2-ratio thresholds that should be used for a sample to call CNA states. Copy neutral is 0. Kernel type is gaussian}
\value{A dataframe with 2 thresholds.}
\author{Emilie Lalonde and Veronica Sabelnykova}
\examples{
get.sample.specific.cna.thresholds(cna.data = rnorm(1000, 5, 4), percent = 0.9);
}
