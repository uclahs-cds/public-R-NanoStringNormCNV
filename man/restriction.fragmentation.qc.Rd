\name{restriction.fragmentation.qc}
\alias{restriction.fragmentation.qc}
\title{Perform QAQC using restriction fragmentation controls}
\description{Perform QAQC using restriction fragmentation controls and create diagnostic figures}
\usage{
	restriction.fragmentation.qc(raw.data)
}
\arguments{
  \item{raw.data}{A gene by sample data-frame of raw NanoString counts of restriction fragmentation probes. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
}
\details{
	Restriction fragmentation controls check for proper denaturation and digestion of DNA in AluI digested samples. The following results are expected: a) probes C and D with counts of > 200 if DNA had correctly denatured prior to digestion and b) probe ratios (C+D / A+B) of > 10 if DNA digestion was complete. The user is warned if low count or low ratio samples are detected.
	\cr Diagnostic plot is saved as a tiff file to the working directory.
	}
\value{A data-frame containing 1) C and D probe mean counts, 2) A and B probe mean counts, and 3) C+D / A+B ratios}
\references{
	See NanoString website for PDFs on analysis guidelines:
		\url{https://www.nanostring.com/support/product-support/support-documentation}

	The NanoString assay is described in the paper:
		Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\author{Cindy Yao and Emilie Lalonde}
\examples{\dontrun{
	# load raw data
	data(NanoString.DNA.raw);

	# run QC
	qc.results <- restriction.fragmentation.qc(raw.data = NanoString.DNA.raw);
}}
