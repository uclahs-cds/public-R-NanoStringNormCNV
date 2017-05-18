\name{positive.control.qc}
\alias{positive.control.qc}
\title{Perform QAQC using positive controls}
\details{Determine the correlation between expected and observed positive control probe values}
\usage{
	positive.control.qc(raw.data)
}
\arguments{
  \item{raw.data}{A gene by sample data-frame of raw NanoString counts for positive control probes. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
}
\description{
	Positive control probes are compared to the expected positive control standard curve, as provided in the NanoString manual, y = 171.23x + 214.12. An R squared value of > 0.95 is expected for every sample and the user is warned if this is not satisfied.
}
\value{A data-frame of R squared values by sample ID}
\references{
	See NanoString website for PDFs on analysis guidelines:
		\url{https://www.nanostring.com/support/product-support/support-documentation}

	The NanoString assay is described in the paper:
		Fortina, P. & Surrey, S. Digital mRNA profiling. Nature Biotechnology 26, 293-294 (2008).
}
\author{Cindy Yao and Emilie Lalonde}
\seealso{\code{\link{make.positive.control.plot}}}
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
