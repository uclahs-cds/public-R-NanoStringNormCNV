\name{call.cna.states}
\alias{call.cna.states}
\title{Convert log2ratio of CNA segments into discrete CNA calls}
\description{Will convert log2ratio CNA calls into CNA states for a cohort of patients}
\usage{call.cna.states(patient.data, cna.thresholds)}
\arguments{
	\item{cna.data}{A data frame/matrix containing the CNA segment information per sample. Must contain a column named log2ratio}
	\item{thresholds}{A 4-item vector specifying homoz loss, het loss, single-copy gain and multi-level amplification, respectively. Order is important!)}
}
\details{Will convert log2-ratio CNA calls into CNA states for a set of samples.}
\value{A vector of copy number states.}
\author{Emilie Lalonde}
