\name{get.cna.thresholds}
\alias{get.cna.thresholds}
\title{Determine thresholds for CNA calling}
\description{Determines tumour/normal ratio thresholds used for calling CNA states. Basic formula is median - x*SD.}
\usage{get.cna.thresholds(ratios, percent = 0.8)}
\arguments{
  \item{ratios}{A vector of tumour/normal ratios.}
  \item{percent}{The percent of density you would like to have between the threshold bounds (0.1 - 0.99). Default is 0.8}
}
\details{This script will return the ratio thresholds that should be used to call CNA states for a given set of tumour/normal ratios. Thresholds are based on the standard deviations of the central 50\% of probes. Kernel type is gaussian.}
\value{A data-frame of 2 thresholds.}
\author{Emilie Lalonde and Veronica Sabelnykova}
\examples{
get.cna.thresholds(ratios = rnorm(1000, 5, 4), percent = 0.9);
}
