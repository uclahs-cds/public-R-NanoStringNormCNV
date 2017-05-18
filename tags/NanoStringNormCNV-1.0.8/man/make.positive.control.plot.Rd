\name{make.positive.control.plot}
\alias{make.positive.control.plot}
\title{Make heatmaps of positive control correlations}
\description{Create heatmaps of correlations between target concentrations and actual counts of positive control probes.}
\usage{
make.positive.control.plot(correlations, covs = NULL, print.x.labels = TRUE)
}
\arguments{
  \item{correlations}{A data-frame of R squared values by sample ID. See output of \code{\link{positive.control.qc}}}
  \item{covs}{A data-frame of sample annotation used to create plot covariates. Recommended to use output of \code{\link{load.phenodata}}. Accepts sample 'Cartridge' and 'Type' information; requires 'SampleID' column to match to sample IDs in 'correlations' variable. Defaults to NULL}
  \item{print.x.labels}{Whether to display the sample IDs on the heatmap's x-axis. Defaults to TRUE}
}
\details{Creates two heatmaps that display the sample R squared values, ascending.}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
\examples{\dontrun{
	#load raw data
	data(NanoString);
	data(PhenoData);

	# run QC
	r.squared <- positive.control.qc(raw.data = NanoString);

	# plot results
	make.positive.control.plot(
		correlations = r.squared,
		covs = PhenoData
		);
}}
