\name{NanoString.norm}
\alias{NanoString.norm}
\docType{data}
\title{Dataset of normalized NanoString counts from prostate cancer patients}
\description{Normalized NanoString counts from 349 probes across 65 NanoString samples. The NanoString samples are the same as those found in the \code{\link{NanoString}} dataset. Invariant probe normalization was applied to raw counts using \code{\link{normalize.global}}. For additional sample information, refer to the \code{\link{PhenoData}} dataset.}
\format{
	A data-frame with 68 columns and 349 rows. The first three columns provide probe information:
	\describe{
		\item{CodeClass:}{probe type which can include: Endogenous, Housekeeping, Invariant, Negative, Positive, RestrictionSite. In this dataset, housekeeping probes were set manually. Probes were either selected from a list of 'Endogenous' probes found to consistently have minimal variance across data or, alternatively, new ones were simulated by adding small noise to existing housekeeping probes (the latter can be identified by the 'SIM' name prefix).}
		\item{Name:}{probe name which is unique to the row}
		\item{Accession:}{probes with matching Accession values correspond to the same gene/segment (for endogenous probes)}
	}
	The following columns indicate the sample (column names are sample IDs). The rows indicate the probe. The contents of the data-frame are raw (unnormalized) NanoString counts.}
\usage{NanoString.norm}
\author{Dorota Sendorek}
\examples{
	data(NanoString.norm);
}
\keyword{datasets}
