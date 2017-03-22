\name{evaluate.replicates}
\alias{evaluate.replicates}
\title{Perform evaluation of replicates}
\description{Evaluate replicate variance and concordance per sample post-CNA calling}
\usage{
evaluate.replicates(phenodata, normalized.data = NULL, cna.rounded = NULL)
}
\arguments{
  \item{phenodata}{A data-frame containing phenotype information. Must include columns 'SampleID', 'Name', 'Type', and 'HasReplicate'. Recommended to use output of \code{\link{load.phenodata}}}
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. Defaults to NULL.}
  \item{cna.rounded}{A gene by sample data-frame of CNA calls. See 'cna.rounded' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL.}
}
\details{A wrapper function for \code{\link{calculate.replicate.variance}} and \code{\link{calculate.replicate.concordance}}. Must provide normalized counts to calculate replicate variance. Must provide CNA calls to calculate replicate concordance. Replicates are identified by their matching sample names (phenodata column 'Name').}
\value{A list of 7 items: 1) normalized count variance scores 2) CNA concordance scores 3) concordance summary for each sample 4) replicate normalized count data 5) phenotype info used to calculate normalized count variance 6) replicate CNA data 7) phenotype info used to calculate CNA concordance. Items are set to NULL if not available due to missing data.}
\author{Emilie Lalonde and Dorota Sendorek}
