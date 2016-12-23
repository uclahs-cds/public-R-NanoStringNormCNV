\name{get.tumour.normal.ratio}
\alias{get.tumour.normal.ratio}
\title{Calculate tumour/normal ratios}
\description{Calculates tumour ratios from sample counts}
\usage{get.tumour.normal.ratio(normalized.data, reference, chip.info, per.chip = FALSE)}
\arguments{
    \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts for tumour and reference samples}
    \item{reference}{A vector of the reference sample name(s)}
    \item{chip.info}{Data-frame of sample annotation with columns 'Cartridge' and 'SampleID'}
    \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE}
}
\details{None}
\value{A gene by sample matrix, where column names are sample IDs and row names are probe names}
\author{Emilie Lalonde}
    

