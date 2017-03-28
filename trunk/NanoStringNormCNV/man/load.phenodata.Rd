\name{load.phenodata}
\alias{load.phenodata}
\title{Load sample annotation and phenotype information}
\description{Load and check the annotation and phenotype information for all tumour and reference samples in the experiment}
\usage{load.phenodata(fname, separator = 'comma')}
\arguments{
  \item{fname}{The filepath for the sample annotation and phenotype information file. Expects a table with the following column headers: 'SampleID' (must be unique to each row/replicate), 'Name' (specific to each tumour or reference sample that has one or more replicates in the study), 'Patient' (specific to each patient that provided one or more tissue samples for the study), 'Type' (options are 'Tumour' or 'Reference'), 'Cartridge' (numeric values only), 'HasReplicate' (whether sample has replicates (1) or not (0)), 'ReferenceID' (sample ID of matched normal; use 'missing' if there is no matched normal, use NA for normal samples). Optionally, user can provide a 'Sex' column (containing values 'M' (for male), 'F' (for female), or NA (if unknown)). If this column is not provided, sample sex is automatically set to NA. Please note that sex information is required to call CNAs on chromosomes X and Y.}
  \item{separator}{The input file's field separator character. Accepts 'comma' (default) or 'tab'.}
}
\details{
	Reads in file and performs a thorough check of information provided to ensure it fits the requirements for downstream analysis functions. 
	\cr Be aware that sample IDs and names must contain only alphanumeric characters. Any non-alphanumeric characters will be replaced with a '.'. Please ensure any modified IDs and names are matched in the raw NanoString data to be used in downstream analysis.
	\cr See the \code{\link{PhenoData}} dataset for an example of expected output. 
	}
\value{A data-frame including the following column headers: 'SampleID', 'Name', 'Patient', 'Cartridge', 'Type', 'ReferenceID', 'HasReplicate', 'Sex'}
\author{Dorota Sendorek and Emilie Lalonde}
