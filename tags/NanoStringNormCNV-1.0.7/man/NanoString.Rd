\name{NanoString}
\alias{NanoString}
\docType{data}
\title{Dataset of raw NanoString counts from prostate cancer patients}
\description{Raw NanoString counts from 349 probes across 65 NanoString samples. The NanoString samples are derived from 46 tumour and normal tissue samples obtained from 31 individual prostate cancer patients. There are matched tumour-normal samples for 12 patients. NanoString replicates exist for 17 tumour/normal tissue samples. For additional sample information, refer to the \code{\link{PhenoData}} dataset.}
\format{
	A data-frame with 68 columns and 349 rows. The first three columns provide probe information:
	\describe{
		\item{CodeClass:}{probe type which can include: Endogenous, Housekeeping, Invariant, Negative, Positive, RestrictionSite. In this dataset, housekeeping probes were set manually. Probes were either selected from a list of 'Endogenous' probes found to consistently have minimal variance across data or, alternatively, new ones were simulated by adding small noise to existing housekeeping probes (the latter can be identified by the 'SIM' name prefix).}
		\item{Name:}{probe name which is unique to the row}
		\item{Accession:}{probes with matching Accession values correspond to the same gene/segment (for endogenous probes)}
	}
	The following columns indicate the sample (column names are sample IDs). The rows indicate the probe. The contents of the data-frame are raw (unnormalized) NanoString counts.}
\usage{NanoString}
\author{Dorota Sendorek}
\examples{
	data(NanoString);
}
\keyword{datasets}
