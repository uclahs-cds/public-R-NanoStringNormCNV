\name{get.tumour.normal.ratio}
\alias{get.tumour.normal.ratio}
\title{Calculate tumour/normal ratios}
\description{Calculates tumour ratios from sample NanoString count data}
\usage{get.tumour.normal.ratio(normalized.data, reference, chip.info, per.chip = FALSE)}
\arguments{
    \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts for tumour and reference samples}
    \item{reference}{A vector of the reference sample name(s)}
    \item{chip.info}{A data-frame of chip information with columns 'Cartridge' and 'SampleID'. See output of \code{\link{load.phenodata}}}
    \item{per.chip}{Whether to use the reference sample(s) specific to a given chip only. Defaults to FALSE}
}
\value{A gene by sample data-frame of tumour/normal ratios}
\author{Emilie Lalonde}
    

