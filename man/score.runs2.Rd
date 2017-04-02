\name{score.runs2}
\alias{score.runs2}
\title{Scoring NanoStringNormCNV runs alternative}
\description{Scoring NanoStringNormCNV runs by fragmentation method}
\usage{score.runs2(normalized.data, cna.rounded, phenodata, cna.normals = NULL)}
\arguments{
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts for both test and control samples. First three columns must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. See output of \code{\link{normalize.global}} or \code{\link{normalize.per.chip}}.}
  \item{cna.rounded}{A gene by sample data-frame of CNA calls for both test and control samples. See 'cna.rounded' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}}
  \item{phenodata}{A data-frame of sample annotation. Must include columns 'SampleID', 'Patient', 'Cartridge' and 'Type'. See output of \code{\link{load.phenodata}}}
  \item{cna.normals}{A gene by sample data-frame of CNA calls for normal samples only. See 'cna.normals' in output of \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL}
}
\details{Produces summary statistics of the given NanoStringNormCNV run}
\value{A list of two sublists: one for 'Alu1' sample scores and the second for 'Sonicate' sample scores. Each sublist holds ARIs for sample chip info and sample type from both normalized count data and copy number calls.}
\author{Dorota Sendorek}
