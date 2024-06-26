\name{generate.plot.legend}
\alias{generate.plot.legend}
\title{Generate plot legend object}
\description{Generate a legend object for plotting purposes}
\usage{
generate.plot.legend(cov.info)
}
\arguments{
  \item{cov.info}{A list of covariate information where each item is a vector of all the sample values for a given covariate. Currently accepts only the following covariates: 'Type', 'Cartridge', 'CodeClass'}
}
\details{Called by \code{\link{visualize.results}}, this function creates a list defining the covariate legend that can be passed onto a plot}
\value{A list where each item contains colour, label and title information for a given covariate}
\author{Dorota Sendorek and Emilie Lalonde}
\seealso{NanoStringNormCNV:::generate.plot.covariates}
\keyword{internal}
