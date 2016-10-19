\name{generate.plot.covariates}
\alias{generate.plot.covariates}
\title{Generate plot covariate object}
\description{Generate a covariate object for plotting purposes}
\usage{
generate.plot.covariates(cov.info)
}
\arguments{
  \item{cov.info}{A sample by covariate data-frame. Currently accepts only the following covariate columns: 'type', 'cartridge', 'CodeClass'}
}
\details{Creates a covariate annotate (as a fully formed list) that can be passed onto a plot}
\value{A list where each item contains border colour and fill colour information for a given covariate}
\author{Dorota Sendorek and Emilie Lalonde}
\seealso{\code{\link{generate.plot.legend}}}
