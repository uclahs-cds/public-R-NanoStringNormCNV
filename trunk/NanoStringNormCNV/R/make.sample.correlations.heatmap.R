
make.sample.correlations.heatmap <- function(nano.counts, cor.method = 'pearson', fname.stem = NULL, covs = NULL) {
	# check and set up covariates
	if (!is.null(covs)) {
		# check completeness and order
		if (! all( colnames(nano.counts) %in% covs$SampleID )) {
			stop("Must provide covariate information for every sample!");
		} else {
			covs <- covs[match(colnames(nano.counts), covs$SampleID),];
			}

		covs <- covs[, !(names(covs) == 'SampleID'), drop = FALSE];
		rownames(covs) <- NULL;

		# create covariate object
		cov.obj <- NanoStringNormCNV::generate.plot.covariates(cov.info = covs);
	} else {
		cov.obj <- NULL;
		}

	# set up legend
	if (!is.null(covs)) {
		covs.legend <- NanoStringNormCNV::generate.plot.legend(cov.info = as.list(covs));
	} else {
		covs.legend <- NULL;
		}

	# set up file name
	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# get the inter-array correlations
	correlations <- cor(x = nano.counts, use = 'all.obs', method = cor.method);

	# plot
	cov.border <- list(col = 'black', lwd = 1.5);
	BoutrosLab.plotting.general::create.heatmap(
		x = correlations,
		filename = paste0(Sys.Date(), fname.stem, '_inter-sample-correlation_heatmap.tiff'),
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
