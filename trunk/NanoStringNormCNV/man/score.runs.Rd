\name{score.runs}
\alias{score.runs}
\title{Scoring runs}
\description{Scoring normalization run using replicate information}
\usage{score.runs(replicates, normalized, cnas, sample.annot, normals = NULL)}
\arguments{
  \item{replicates}{A list of replicate information. See output in \code{\link{evaluate.replicates}}}
  \item{normalized}{A gene by sample data-frame of normalized NanoString counts}
  \item{cnas}{A gene by sample data-frame of CNA calls for test and control samples. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs. See 'cna.rounded' in output of \code{\link{call.cnas.with.matched.normals}} or \code{\link{call.cnas.with.pooled.normals}}}
  \item{sample.annot}{A data-frame containing phenotype information. Must include columns 'SampleID', 'Patient', 'Cartridge' and 'Type'. See output in \code{\link{load.phenodata}}}
  \item{normals}{A gene by sample data-frame of CNA calls for normal samples only. See 'cna.normals' in output of \code{\link{call.cnas.with.pooled.normals}}. Defaults to NULL}
}
\details{Produces summary statistics of the given NanoStringNormCNV run}
\value{A list of the following items: 1) ARI for replicate all sample count correlations, replicate tumour only count correlations, tissue type, and CNA calls 2) standard deviation for invariant probes, housekeeping probes, and both 3) mean replicate concordance 4) proportion of discordant genes 5) number of normals with CNAs 6) total of CNAs in normals 7) total of CNAs overall }
\author{Emilie Lalonde}
