\name{get.ari}
\alias{get.ari}
\title{Determine the adjusted rand index of a given feature}
\description{Determine the adjusted rand index of a given feature for continuous or discrete data}
\usage{
get.ari(data, feature, discrete.data = TRUE)
}
\arguments{
  \item{data}{The data-frame to be clustered}
  \item{feature}{A vector with the variable of interest. Please note, it will be converted into a factor in the function.}
  \item{discrete.data}{A boolean indicating whether data is discrete. If TRUE (default), Jaccard distance and ward clustering will be used and, if FALSE, Pearson correlation and complete clustering is used.}
}
\details{Determine the adjusted rand index of a variable given continuous or discrete data to be clustered.}
\value{The adjusted rand index of feature}
\author{Cindy Yao and Emilie Lalonde}
