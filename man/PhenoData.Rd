\name{PhenoData}
\alias{PhenoData}
\docType{data}
\title{Dataset of sample annotation and phenotype information for 65 prostate cancer samples}
\description{Sample annotation and phenotype information describing 65 NanoString prostate cancer samples. The same NanoString samples are found in the \code{\link{NanoString}} dataset.}
\format{
	A data-frame with 9 columns and 65 rows. Each row indicates a different NanoString sample. The columns either annotate the sample or describe a feature of it:
	\describe{
		\item{SampleID}{The sample ID is unique to each row (or NanoString sample).}
		\item{Patient}{The patient ID is specific to an individual patient that provided one or more tissue samples for the study.}
		\item{Name}{The sample name is unique to a tumour/normal tissue sample. Each tissue sample can have one or more replicates (also referred to as NanoString samples) which are distinguished by sample IDs ('SampleID').}
		\item{Cartridge}{The cartridge the NanoString sample was processed on. All numeric values.}
		\item{Type}{The sample type, either 'Tumour' or 'Reference' (also referred to as a normal).}
		\item{ReferenceID}{The sample ID of the sample's matched normal. This only applies to tumour samples so for reference samples this value is NA. For tumour samples with no matched normal, 'missing' is used.}
		\item{HasReplicate}{Whether the sample has replicates (1) or not (0).}
		\item{Sex}{The sex of the patient that provided the sample. The following notation is used: 'M' for male, 'F' for female, NA for unknown. Sample sex is required to perform CNA calling on sex chromosome probes.}
		\item{Fragmentation}{The DNA fragmentation used: sonication or AluI restriction site digestion. This column is optional. The information is only used when selecting samples for a quality check using restriction fragmentation controls.}
	}
}
\usage{PhenoData}
\author{Dorota Sendorek}
\examples{
	data(PhenoData);
}
\keyword{datasets}
