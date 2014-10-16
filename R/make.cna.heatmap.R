
make.cna.heatmap <- function(nano.cnas, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL, rounded = FALSE, clust.dim = 'both'){

	dist.method <- ifelse(rounded, 'jaccard', 'correlation');
	create.heatmap(
		x = nano.cnas,
		filename = generate.filename(fname.stem, 'cna_heatmap', 'png'),
		cluster.dimensions = clust.dim,
		clustering.method = 'centroid',
		xlab.label = 'Genes',
		ylab.cex = 2,
		ylab.label = 'Samples',
		xlab.cex = 2,
		rows.distance.method = dist.method,
		cols.distance.method = dist.method,
		colour.centering.value = 3,
		at = seq(1,max(nano.cnas),by = 0.1),
		colourkey.labels = c(0,2,4),
		colourkey.labels.at = c(1,3,5),
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = c('blue', 'white', 'red'),
		resolution = 600
		);

	}
