\name{NanoString.DNA.norm}
\alias{NanoString.DNA.norm}
\docType{data}
\title{Dataset of normalized NanoString DNA counts from prostate cancer patients}
\description{
Normalized NanoString DNA counts from 349 probes across 65 NanoString samples. The NanoString samples are the same as those found in the \code{\link{NanoString}} dataset. Invariant probe normalization was applied to raw counts using \code{\link{normalize.global}}. For additional sample information, refer to the \code{\link{PhenoData}} dataset.
}
\format{
	A data-frame with 68 columns and 349 rows. The first three columns provide probe information:
	\describe{
		\item{CodeClass:}{probe type which can include: Endogenous, Housekeeping, Invariant, Negative, Positive, RestrictionSite.}
		\item{Name:}{probe name which is unique to the row}
		\item{Accession:}{probes with matching Accession values correspond to the same gene/segment (for endogenous probes)}
	}
	The following columns indicate the sample (column names are sample IDs). The rows indicate the probe. The contents of the data-frame are normalized NanoString DNA counts.

	Sample ID notation details are as follows. The first 8 characters are the patient ID. Sample type is denoted by one of the following: F = tumour (fresh frozen), P = tumour (FFPE) or B = reference (blood). Additionaly, 'M1' and 'M2' are used here to differentiate between replicates.
}
\usage{NanoString.DNA.norm}
\author{Dorota Sendorek}
\examples{
	data(NanoString.DNA.norm);
}
\keyword{datasets}
