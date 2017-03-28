\name{generate.plot.covariates}
\alias{generate.plot.covariates}
\title{Generate plot covariate object}
\description{Generate a covariate object for plotting purposes}
\usage{
generate.plot.covariates(plotting.data, sample.covariates = NULL, gene.covariates = NULL)
}
\arguments{
  \item{plotting.data}{A matrix or data-frame that contains the main plotting data. Sample IDs are expected for column names}
  \item{sample.covariates}{A sample by covariate data-frame. Currently only accepts sample covariates 'Type' and 'Cartridge'. Must also contain a 'SampleID' column to match to 'nano.counts' sample IDs. Defaults to NULL}
  \item{gene.covariates}{A gene by covariate data-frame. Currently only accepts gene covariate 'CodeClass'. A 'Name' column is mandatory to match to 'nano.counts' gene probe names. Defaults to NULL}
}
\details{Called by \code{\link{visualize.results}}, this function creates a covariate annotate (as a fully formed list) that can be passed onto a plot}
\value{A list of 2 sub-lists that contain border colour and fill colour information for 1) sample covariates and 2) gene covariates}
\author{Dorota Sendorek and Emilie Lalonde}
\seealso{\code{\link{generate.plot.legend}}}
