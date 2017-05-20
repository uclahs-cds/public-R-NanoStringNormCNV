\name{make.counts.heatmap}
\alias{make.counts.heatmap}
\title{Make a heatmap for NanoString counts}
\description{Make a heatmap of counts}
\usage{
make.counts.heatmap(
  nano.counts, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL,
  clust.dim = 'both', clust.method = 'euclidean', print.ylab = NULL
  )
}
\arguments{
  \item{nano.counts}{A gene by sample matrix or data-frame of NanoString counts, where row names are probe names and column names are sample IDs. Counts can be normalized or raw.}
  \item{fname.stem}{To use in filename. Description of how CNAs were called is suggested. Defaults to NULL}
  \item{covs.rows}{A sample by covariate data-frame. Currently only accepts sample covariates 'Type' and 'Cartridge' (see output of \code{\link{load.phenodata}}). A 'SampleID' column is mandatory to match to 'nano.counts' sample IDs. Defaults to NULL}
  \item{covs.cols}{A gene by covariate data-frame. Currently only accepts gene covariate 'CodeClass'. A 'Name' column is mandatory to match to 'nano.counts' gene probe names. Defaults to NULL}
  \item{clust.dim}{Which dimensions to cluster. Accepts 'none', 'rows', 'columns', or 'both' (default)}
  \item{clust.method}{Which distance method should be used for clustering. Accepts 'euclidean' (default), 'correlation', 'jaccard', and other methods supported in ?dist}
  \item{print.ylab}{Whether to print the yaxis labels. Passed directly to the BoutrosLab.plotting.general::create.heatmap function. Leave as NULL (default) if you don't want them printed}
}
\details{Make a clustered heatmap of counts for all samples and genes}
\value{None}
\author{Cindy Yao and Emilie Lalonde}
\examples{\dontrun{
  # load data
  data(NanoString);
  data(NanoString.norm);
  data(PhenoData);

  # plot normalized NanoString counts
  make.counts.heatmap(
    nano.counts = NanoString.norm[, -(1:3)],
    fname.stem = 'normalized',
    covs.rows = PhenoData[, c('SampleID', 'Type', 'Cartridge')],
    covs.cols = NanoString[, c('Name', 'CodeClass')]
    ); 

  # plot raw NanoString counts
  # make sure raw count data frame has gene names for rownames!
  NanoString.formatted <- NanoString[, -(1:3)];
  rownames(NanoString.formatted) <- NanoString$Name;

  make.counts.heatmap(
    nano.counts = NanoString.formatted,
    fname.stem = 'raw',
    covs.rows = PhenoData[, c('SampleID', 'Type', 'Cartridge')],
    covs.cols = NanoString[, c('Name', 'CodeClass')]
    );
}}