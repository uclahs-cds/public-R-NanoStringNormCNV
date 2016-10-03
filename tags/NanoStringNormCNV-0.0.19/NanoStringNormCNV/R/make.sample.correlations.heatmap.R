
make.sample.correlations.heatmap <- function(data, cor.method = 'pearson', fname.stem = NULL, covs = NULL, covs.legend = NULL) {
	# get the inter-array correlations and plot
	correlations <- cor(x = data, use = 'all.obs', method = cor.method);

	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	BoutrosLab.plotting.general::create.heatmap(
		x = correlations,
		filename = paste0(Sys.Date(), fname.stem, '_inter-sample-correlation_heatmap.tiff'),
		clustering.method = 'complete',
		plot.dendrograms = 'both',
		xaxis.rot = 90,
		covariates = covs,
		covariates.top = covs,
		covariate.legends = covs.legend,
		legend.title.just = 'left',
		legend.title.fontface = 'bold',
		scale.data = FALSE,
		colourkey.cex = 2.0,
		resolution = 600,
		width = 7
		);
	}
