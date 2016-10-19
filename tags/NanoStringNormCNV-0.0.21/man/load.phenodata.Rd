\name{load.phenodata}
\alias{load.phenodata}
\title{Load phenotype information}
\description{Load and check the phenotype information for NanoString samples}
\usage{load.phenodata(fname)}
\arguments{
  \item{fname}{The filename for phenotype information file. Expects a CSV file with the following column headers: SampleID, Name, Patient, type ('Tumour' or 'Reference'), cartridge (number), has.repl (whether sample has replicate (1) or not (0)), ref.name (name of matched normal, 'missing' if matched normal does not exist, NA for normal samples). Optionally, user can provide a 'sex' column (containing values 'M' (male), 'F' (female), or NA). Please note that sex information is required for calling CNAs on chromosomes X/Y.}
}
\details{Reads in phenotype information and performs a thorough check of information provided to ensure it fits requirements of downstream analysis functions.
Be aware that sample IDs and names that do not follow the expected naming convention (only alphanumeric characters or '.') will be modified. Please ensure IDs and names match NanoString samples.}
\value{A sample by annotation data-frame, containing the following column headers: SampleID, Name, Patient, cartridge, type, ref.name, has.repl, Chip, sex}
\author{Dorota Sendorek and Emilie Lalonde}
