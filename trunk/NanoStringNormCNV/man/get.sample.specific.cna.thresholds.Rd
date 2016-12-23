\name{get.sample.specific.cna.thresholds}
\alias{get.sample.specific.cna.thresholds}
\title{Determine ratio thresholds for CNA calling}
\description{Determines sample-specific thresholds to call CNA states. Basic formula is median - x*SD.}
\usage{get.sample.specific.cna.thresholds(ratios, percent = 0.8)}
\arguments{
  \item{ratios}{A vector of tumour/normal ratios.}
  \item{percent}{The percent of density you would like to have between the threshold bounds (0.1 - 0.99). Default is 0.8}
}
\details{This script will return the ratio thresholds that should be used to call CNA states in a given sample. Thresholds are based on the standard deviations of the central 50\% of probes. Copy neutral is 2. Kernel type is gaussian.}
\value{A data-frame with 2 thresholds.}
\author{Emilie Lalonde and Veronica Sabelnykova}
\examples{
get.sample.specific.cna.thresholds(ratios = rnorm(1000, 5, 4), percent = 0.9);
}
