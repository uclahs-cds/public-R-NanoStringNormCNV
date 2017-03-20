\name{invariant.probe.norm}
\alias{invariant.probe.norm}
\title{Normalize data using invariant probes}
\description{Normalize data using invariant probes to control for differences in sample input}
\usage{
invariant.probe.norm(nanostring.data, phenodata = NULL)
}
\arguments{
  \item{nanostring.data}{A data-frame of NanoString counts that are either raw or have already been normalized by a different normalization technique. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{A data-frame with fields 'SampleID' and 'Type' (specifying whether 'Tumour' or 'Reference'). Defaults to NULL}
}
\details{
	NanoString counts are multiplied by a normalization factor (mean INV probe counts in sample / mean(mean INV probe count of each sample). This procedure is not carried out on 'Positive', 'Negative' and 'RestrictionSite' code class probes. Please note that invariant probe counts that are < 100 signify low DNA input. This may lead to unreliable CNA calls, especially if the low counts are in the reference samples! 
	This function is called from \code{\link{normalize.global}} and \code{\link{normalize.per.chip}}, if invariant probe normalization is requested.
}
\value{A gene by sample data-frame of normalized counts, where first three columns contain 'CodeClass', 'Name', and 'Accession' information}
\references{
	See NanoString website for PDFs on analysis guidelines:
		\url{http://www.nanostring.com/support/prod-lit/}

	The NanoString assay is described in the paper:
		Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::sample.content.normalization}
