\name{get.tumour.normal.ratio}
\alias{get.tumour.normal.ratio}
\title{Calculate tumour/normal ratios}
\description{Calculates tumour ratios from sample counts}
\usage{get.tumour.normal.ratio(ns.counts, ref, chips.info, per.chip = FALSE)}
\arguments{
    \item{ns.counts}{A gene by sample data-frame of raw NanoString counts for tumour and reference samples}
    \item{ref}{A vector of the reference sample name(s)}
    \item{chips.info}{Data-frame of sample annotation with columns 'Cartridge' and 'SampleID'}
    \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE}
}
\details{None}
\value{A gene by sample matrix, where column names are sample IDs and row names are probe names}
\author{Emilie Lalonde}
    

