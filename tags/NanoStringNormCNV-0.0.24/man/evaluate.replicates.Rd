\name{evaluate.replicates}
\alias{evaluate.replicates}
\title{Perform evaluation of replicates}
\description{Evaluate replicate variance and concordance per sample}
\usage{
evaluate.replicates(normalized.data, phenodata, cnas)
}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts}
  \item{phenodata}{Data-frame containing phenotype information (specifically 'SampleID', 'Name', 'type', and 'has.repl'). See output in \code{\link{load.phenodata}}}
  \item{cnas}{A gene by sample data-frame of NanoString CNA calls, where column names are sample IDs and row names are gene probe names}
}
\details{A wrapper function for \code{\link{calculate.replicate.variance}} and \code{\link{calculate.replicate.concordance}}}
\value{A list of 7 items: 1) normalized count variance scores 2) CNA call concordance scores 3) concordance summary for each sample name 4) replicate normalized count data 5) phenotype info used to calculate normalized count variance 6) replicate NanoString CNA calls 7) phenotype info used to calculate CNA call concordance}
\author{Emilie Lalonde and Dorota Sendorek}
