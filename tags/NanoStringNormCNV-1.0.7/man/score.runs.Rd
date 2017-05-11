\name{score.runs}
\alias{score.runs}
\title{Scoring NanoStringNormCNV runs}
\description{Scoring NanoStringNormCNV runs using replicate information}
\usage{score.runs(replicate.eval, normalized.data, cna.rounded, phenodata, cna.normals = NULL)}
\arguments{
  \item{replicate.eval}{A list of replicate evaluation information. Use output of \code{\link{evaluate.replicates}}}
  \item{normalized.data}{A gene by sample data-frame of normalized NanoString counts for both test and control samples. First three columns must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. See output of \code{\link{normalize.global}} or \code{\link{normalize.per.chip}}.}
  \item{cna.rounded}{A gene by sample data-frame of CNA calls for both test and control samples. See 'cna.rounded' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}}
  \item{phenodata}{A data-frame of sample annotation. Must include columns 'SampleID', 'Patient', 'Cartridge' and 'Type'. See output of \code{\link{load.phenodata}}}
  \item{cna.normals}{A gene by sample data-frame of CNA calls for normal samples only. See 'cna.normals' in output of \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL}
}
\details{Produces summary statistics of the given NanoStringNormCNV run}
\value{A list of two lists. The first sublist contains scores. Adjusted Rand indices (ARI) for the following features: [1] cartridge, [2] patients (using normalized counts), [3] tissue type, and [4] patients (using CNA calls); normalized count standard deviations for [5] invariant probes, [6] housekeeping probes, and [7] invariant and housekeeping probes combined; [8] mean replicate concordance; [9] proportion of discordant genes; [10] number of normal samples with one or more CNAs (if applicable); [11] total number CNAs detected across all normal samples (if applicable); [12] total number of CNAs detected across all tumour samples. The second sublist contains the numbers of samples used to calculate the scores.}
\author{Emilie Lalonde}
