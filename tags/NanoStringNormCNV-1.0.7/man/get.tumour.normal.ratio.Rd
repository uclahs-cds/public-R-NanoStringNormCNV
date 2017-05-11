\name{get.tumour.normal.ratio}
\alias{get.tumour.normal.ratio}
\title{Calculate tumour/normal ratios}
\description{Calculates tumour/normal ratios using NanoString count data}
\usage{get.tumour.normal.ratio(normalized.data, reference, chip.info = NULL, per.chip = FALSE)}
\arguments{
    \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts of tumour and reference samples. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs.}
    \item{reference}{A vector of the reference sample name(s)}
    \item{chip.info}{A data-frame of chip information for each sample in 'normalized.data', including columns 'Cartridge' and 'SampleID', or NULL (default). Required if 'per.chip' = TRUE. Recommended to use output of \code{\link{load.phenodata}}}
    \item{per.chip}{Whether to use only the reference sample(s) specific to a given chip. Defaults to FALSE.}
}
\value{A gene by sample data-frame of tumour/normal ratios}
\author{Emilie Lalonde}
    

