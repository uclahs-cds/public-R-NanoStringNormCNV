\name{invariant.probe.norm}
\alias{invariant.probe.norm}
\title{Normalize data using invariant probes}
\description{Normalize data using invariant (INV) probes to control for differences in sample input}
\usage{
invariant.probe.norm(nanostring.data, phenodata = NULL)
}
\arguments{
  \item{nanostring.data}{A gene by sample data-frame of NanoString counts (can either be raw or have already been normalized by a separate normalization technique). First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{A data-frame of sample annotation with fields 'SampleID' and 'Type' (specifying whether sample is 'Tumour' or 'Reference'). Recommended to use output of \code{\link{load.phenodata}}. Defaults to NULL}
}
\details{
	NanoString counts are multiplied by a normalization factor (mean(mean INV probe count of each sample) / mean INV probe count of given sample). This procedure is not carried out on 'Positive', 'Negative' and 'RestrictionSite' code class probes.
	\cr Be aware that invariant probe counts of < 100 signify low DNA input. This may lead to unreliable CNA calls, especially if the low counts are in the reference samples! 
	\cr This function is called from \code{\link{normalize.global}} and \code{\link{normalize.per.chip}}, if invariant probe normalization is requested.
	\cr Diagnostic plots are saved as tiff files to the working directory.
}
\value{A gene by sample data-frame of normalized counts, where the first three columns contain 'CodeClass', 'Name', and 'Accession' information}
\references{
	See NanoString website for PDFs on analysis guidelines:
		\url{https://www.nanostring.com/support/product-support/support-documentation}

	The NanoString assay is described in the paper:
		Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::sample.content.normalization}
\examples{\dontrun{
	# load raw data
	data(NanoString);
	data(PhenoData)

	# normalize
	normalized.data <- invariant.probe.norm(
		nanostring.data = NanoString,
		phenodata = PhenoData
		);	
}}
