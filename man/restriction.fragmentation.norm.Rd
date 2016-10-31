\name{restriction.fragmentation.norm}
\alias{restriction.fragmentation.norm}
\title{Do QAQC based on restriction fragmentation probes}
\description{Do QAQC based on restriction fragmentation probes and make diagnostic figures}
\usage{
restriction.fragmentation.norm(nano.df)
}
\arguments{
  \item{nano.df}{A gene by sample data-frame of raw NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
}
\details{Do QAQC based on restriction fragmentation probes and make diagnostic figures. For AluI digested DNA samples, expect the following: a) probes C and D with counts > 200 and b) probes C+D / A+B > 10}
\value{A data-frame containing 1) C and D probe mean counts, 2) A and B probe mean counts, and 3) C+D / A+B ratios}
\references{NanoString manual p 7}
\author{Cindy Yao and Emilie Lalonde}
