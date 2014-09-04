
make.counts.heatmap <- function(nano.counts, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL){

	create.heatmap(
		x = nano.counts,
		filename = generate.filename(fname.stem, 'counts_heatmap', 'png'),
		xlab.label = 'Genes',
		ylab.cex = 2,
		ylab.label = 'Samples',
		xlab.cex = 2,
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = default.colours(2, palette.type = 'qual'),
		resolution = 600
		);

	}
