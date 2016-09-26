\name{get.ari}
\alias{get.ari}
\title{Determine the adjusted rand index of a given feature}
\description{Determine the adjusted rand index of a given feature for continuous or discrete data}
\usage{
get.ari(data, feature, discrete.data = FALSE)
}
\arguments{
  \item{data}{The data to be clustered}
  \item{feature}{A vector with the variable of interest. Will be converted into a factor in the function.}
  \item{discrete.data}{A boolean indicating whether data is discrete, defaults to FALSE. If TRUE, Jaccard distance and ward clustering will be used instead of Pearson correlation and complete clustering.}
}
\details{Determine the adjusted rand index of a variable given continuous or discrete data to be clustered.}
\value{The adjusted rand index of feature}
\author{Cindy Yao and Emilie Lalonde}
