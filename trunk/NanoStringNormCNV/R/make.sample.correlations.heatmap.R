make.sample.correlations.heatmap <- function(nano.counts, cor.method = 'pearson', fname.stem = NULL, covs = NULL) {
	# set up legend
	if (!is.null(covs)) {
		# get required samples
		covs <- covs[covs$SampleID %in% names(nano.counts),];

		# remove 'Type' if all samples are either 'Tumour' or 'Reference'
		if ("Type" %in% colnames(covs)) {
			if (nlevels(as.factor(covs$Type)) == 1) {
				covs <- covs[, !colnames(covs) %in% "Type", drop = FALSE];
				}
			}

		# covariates
		cov.objs <- NanoStringNormCNV::generate.plot.covariates(
			plotting.data = nano.counts,
			sample.covariates = covs
			);
		cov.obj <- cov.objs[['sample']];

		# legend
		covs <- covs[, names(covs) != 'SampleID', drop = FALSE];
		covs.legend <- NanoStringNormCNV::generate.plot.legend(cov.info = as.list(covs));
	} else {
		covs.legend <- NULL;
		cov.obj <- NULL;
		}

	# set up file name
	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# get the inter-array correlations
	correlations <- cor(x = nano.counts, use = 'all.obs', method = cor.method);

	# plot
	cov.border <- list(col = 'black', lwd = 1.5);
	BoutrosLab.plotting.general::create.heatmap(
		x = correlations,
		filename = paste0(Sys.Date(), fname.stem, '_inter-sample-correlation-heatmap.tiff'),
		clustering.method = 'complete',
		plot.dendrograms = 'both',
		xaxis.rot = 90,
		xaxis.tck = 0,
		yaxis.tck = 0,
		covariates = cov.obj,
		covariates.top.grid.border = cov.border,
		covariates.top = cov.obj,
		covariates.grid.border = cov.border,
		covariate.legends = covs.legend,
		legend.title.just = 'left',
		legend.title.fontface = 'bold',
		scale.data = FALSE,
		colourkey.cex = 2.0,
		resolution = 600,
		width = 7
		);
	}
