\name{invariant.probe.norm}
\alias{invariant.probe.norm}
\title{Normalize data using invariant probes}
\description{Normalize data using invariant probes to control for differences in sample input}
\usage{
invariant.probe.norm(normalized.data, phenodata = NULL)
}
\arguments{
  \item{normalized.data}{A data-frame of normalized NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
  \item{phenodata}{Data-frame with fields 'SampleID' and 'Type' (specifying whether Tumour or Reference). Defaults to NULL}
}
\details{Normalize data using invariant probes to control for differences in sample input. Counts are multiplied by a normalization factor (mean INV probe counts in sample / mean(mean INV probe count of each sample). Invariant probe counts that are < 100 signify low DNA input. This may lead to unreliable CNA calls, especially if the low counts are in the reference samples!}
\value{Nanostring data frame with normalized counts}
\references{Manual pp 8-10}
\author{Cindy Yao and Emilie Lalonde}
\seealso{NanoStringNorm::sample.content.normalization}
