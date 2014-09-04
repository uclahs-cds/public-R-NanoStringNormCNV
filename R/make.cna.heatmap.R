
make.cna.heatmap <- function(nano.cnas, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL){

	create.heatmap(
		x = nano.cnas,
		filename = generate.filename(fname.stem, 'cna_heatmap', 'png'),
		xlab.label = 'Genes',
		ylab.cex = 2,
		ylab.label = 'Samples',
		xlab.cex = 2,
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = c('blue', 'white', 'red'),
		resolution = 600
		);

	}
