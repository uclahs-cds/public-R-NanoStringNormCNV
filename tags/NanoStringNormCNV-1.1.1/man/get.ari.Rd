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
\examples{\dontrun{
	# load data
	data(NanoString.norm);
	data(PhenoData);

	# call CNAs
	cnas <- call.cnas.with.pooled.normals(
		normalized.data = NanoString.norm,
		phenodata = PhenoData
		);

	# evaluate results using replicates
	evaluation <- evaluate.replicates(
		phenodata = PhenoData,
		normalized.data = NanoString.norm,
		cna.rounded = cnas$rounded
		);
	
	## example 1
	# determine how well the copy number calls cluster by sample patient
	patient.ari <- get.ari(
		data.to.cluster = evaluation$cna.calls,
		feature = PhenoData[
			match(
				colnames(evaluation$cna.calls),
				PhenoData$SampleID
				),
			]$Patient,
		is.discrete = TRUE
		);

	## example 2
	# determine to what extent the normalized counts cluster by sample cartridge
	# log values, if appropriate
	if (all(unlist(NanoString.norm) >= 0)) {
	    count.data <- log10(NanoString.norm[, -c(1:3)] + 1);
	} else {
	    count.data <- NanoString.norm[, -c(1:3)];
	    }

	cartridge.ari <- get.ari(
	    data.to.cluster = count.data,
	    feature = PhenoData$Cartridge[
	    	match(
	    		colnames(NanoString.norm[, -(1:3)]),
	    		PhenoData$SampleID
	    		)
	    	],
	    is.discrete = FALSE
	    );
}}
