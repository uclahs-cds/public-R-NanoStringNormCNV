\name{load.phenodata}
\alias{load.phenodata}
\title{Load phenotype information}
\description{Load and check the phenotype information for all samples}
\usage{load.phenodata(fname, separator = 'comma')}
\arguments{
  \item{fname}{The filepath for the phenotype information file. Expects a table with the following column headers: SampleID, Name, Patient, Type ('Tumour' or 'Reference'), Cartridge (number only), HasReplicate (whether sample has replicate (1) or not (0)), ReferenceID (name of matched normal, 'missing' if matched normal does not exist, NA for normal samples). Optionally, user can provide a 'Sex' column (containing values M (male), F (female), or NA). If not provided, sex is automatically set to NA. Please note that sex information is required for calling CNAs on chromosomes X/Y.}
  \item{separator}{The field separator character. Accepts 'comma' (default) or 'tab'.}
}
\details{Reads in phenotype information and performs a thorough check of information provided to ensure it fits requirements of downstream analysis functions. Be aware that sample IDs and names that do not follow the expected naming convention (only alphanumeric characters or '.') will be modified. Please ensure these IDs and names match those in NanoString sample data post-processing.}
\value{A sample by annotation data-frame, containing at least the following column headers: 'SampleID', 'Name', 'Patient', 'Cartridge', 'Type', 'ReferenceID', 'HasReplicate', 'Sex'}
\author{Dorota Sendorek and Emilie Lalonde}
