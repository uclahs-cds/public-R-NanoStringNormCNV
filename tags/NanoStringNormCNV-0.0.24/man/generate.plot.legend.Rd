\name{generate.plot.legend}
\alias{generate.plot.legend}
\title{Generate plot legend object}
\description{Generate a legend object for plotting purposes}
\usage{
generate.plot.legend(cov.info)
}
\arguments{
  \item{cov.info}{A list where each item is a vector containing covariate information. Currently accepts only the following covariates: 'type', 'cartridge', 'CodeClass'}
}
\details{Creates a list defining the covariate legend that can be passed onto a plot}
\value{A list where each item contains colour, label and title information for a given covariate}
\author{Dorota Sendorek and Emilie Lalonde}
\seealso{\code{\link{generate.plot.covariates}}}
