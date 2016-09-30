\name{score.runs}
\alias{score.runs}
\title{Scoring runs}
\description{Scoring normalization run using replicate information}
\usage{score.runs(replicates, normalized, cnas, sample.annot, genes, normals)}
\arguments{
  \item{replicates}{A list of replicate information. Must include items 'norm.counts' and 'raw.counts'}
  \item{normalized}{The NanoString normalized counts data frame}
  \item{cnas}{A data-frame of copy number states for test and control samples}
  \item{sample.annot}{A data-frame with annotation information for each sample. Must include 'cartridge' and 'type' columns}
  \item{genes}{A data-frame of gene probe information. One row per probe containing 'CodeClass', probe name ('Name') and gene name ('Symbol')}
  \item{normals}{A data-frame of copy number states for normal samples only}
}
\details{Produces summary statistics of the given NanoStringNormCNV run}
\value{A list of the following items: ARI for replicate all sample count correlations, replicate tumour only count correlations, sample type, and CNA calls; standard deviation for invariant probes, housekeeping probes, and both; mean replicate concordance; proportion of discordant genes; number of normals with CNAs; total of CNAs in normals; total of CNAs overall }
\author{Emilie Lalonde}
