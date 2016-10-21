\name{get.tumour.normal.ratio}
\alias{get.tumour.normal.ratio}
\title{Calculate tumour/normal ratios}
\description{Calculates tumour ratios from sample counts}
\usage{get.tumour.normal.ratio(ns.counts, ref, chips.info, per.chip = FALSE)}
\arguments{
    \item{ns.counts}{The data associated with the samples}
    \item{ref}{The reference sample(s)}
    \item{chips.info}{Data-frame comprising of the samples' annotations. Required when per.chip = TRUE}
    \item{per.chip}{Whether to use the reference specific to each chip only. Defaults to FALSE}
}
\details{None}
\value{None}
\author{Emilie Lalonde}
    

