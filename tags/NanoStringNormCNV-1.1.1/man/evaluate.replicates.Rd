\name{evaluate.replicates}
\alias{evaluate.replicates}
\title{Perform evaluation using replicates}
\description{Evaluate replicate count variance and CNA call concordance, per sample}
\usage{
evaluate.replicates(phenodata, normalized.data = NULL, cna.rounded = NULL)
}
\arguments{
  \item{phenodata}{A data-frame of sample annotation. Must include columns 'SampleID', 'Name', 'Type', and 'HasReplicate'. Recommended to use output of \code{\link{load.phenodata}}}
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Defaults to NULL}
  \item{cna.rounded}{A gene by sample data-frame of CNA calls. See 'cna.rounded' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL.}
}
\details{Replicates belonging to a single sample are identified as such by matching sample names (phenodata column 'Name'). Normalized counts ('normalized.data') must be provided to calculate replicate variance. CNA calls ('cna.rounded') must be provided to calculate replicate concordance.}
\value{A list of 7 items: 1) normalized count variance scores 2) CNA concordance scores 3) concordance summary for each sample 4) replicate normalized count data used 5) phenotype info used to calculate normalized count variance 6) replicate CNA data used and 7) phenotype info used to calculate CNA concordance. Items are set to NULL if not available due to missing data.}
\author{Emilie Lalonde and Dorota Sendorek}
\examples{\dontrun{
	# load data
	data(NanoString.norm);
	data(PhenoData);

	# call CNAs
	cnas <- call.cnas.with.pooled.normals(
		normalized.data = NanoString.norm,
		phenodata = PhenoData
		);

	# evaluate results using replicates
	evaluation <- evaluate.replicates(
		phenodata = PhenoData,
		normalized.data = NanoString.norm,
		cna.rounded = cnas$rounded
		);
}}
