\name{positive.control.norm}
\alias{positive.control.norm}
\title{Do QAQC based on positive controls}
\details{Find the correlation of positive control probes and the expected amounts}
\usage{
positive.control.norm(nano.df)
}
\arguments{
  \item{nano.df}{A gene by sample data-frame of raw NanoString counts. First three column names must be 'CodeClass', 'Name', and 'Accession', followed by sample IDs}
}
\description{Find the correlation of positive control probes and the expected amounts based on NanoString manual (y = 171.23x + 214.12)}
\value{A data-frame of Pearson correlations by sample}
\author{Cindy Yao and Emilie Lalonde}
