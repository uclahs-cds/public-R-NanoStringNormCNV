
make.samples.correlation.heatmap <- function(data, cor.method = 'pearson', fname, covs, covs.legend){
	# get the inter-array correlations and plot
	correlations <- cor(x = data, use = 'all.obs', method = cor.method);

	BoutrosLab.plotting.general::create.heatmap(
		x = correlations,
		filename = fname,
		clustering.method = 'complete',
		plot.dendrograms = 'both',
		covariates = covs,
		covariates.top = covs,
		covariate.legends = covs.legend,
		legend.title.just = 'left',
		legend.title.fontface = 'bold',
		scale.data = FALSE,
		colourkey.cex = 2.0,
		xaxis.rot = 90,
		resolution = 600,
		width = 7
		);
	}
