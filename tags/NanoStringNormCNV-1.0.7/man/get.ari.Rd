\name{get.ari}
\alias{get.ari}
\title{Determine the adjusted rand index of a given feature}
\description{Determine the adjusted rand index (ARI) of a given feature for continuous or discrete data}
\usage{
get.ari(data.to.cluster, feature, is.discrete = TRUE)
}
\arguments{
  \item{data.to.cluster}{A gene by sample data-frame to be clustered}
  \item{feature}{A vector with the variable of interest. Please note, it will be converted into a factor in the function. If the number of levels (unique elements) is less than 2, ARI will not be calculated.}
  \item{is.discrete}{Whether or not 'data.to.cluster' is discrete. If TRUE (default), Jaccard distance and ward clustering will be used and, if FALSE, Pearson correlation and complete clustering is used.}
}
\details{Determine the adjusted rand index of a feature using clustered data (continuous or discrete).}
\value{The adjusted rand index of feature}
\author{Cindy Yao and Emilie Lalonde}
