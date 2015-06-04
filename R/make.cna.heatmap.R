
make.cna.heatmap <- function(
	nano.cnas, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL,
	rounded = FALSE, clust.dim = 'both', ...
	) {

	# define parameters up front
	dist.method <- ifelse(rounded, 'jaccard', 'euclidean');
	plot.at <- seq(floor(min(nano.cnas, na.rm = TRUE)), ceiling(max(nano.cnas, na.rm = TRUE)), 0.1);
	plot.seq <- seq(min(plot.at), max(plot.at));

	create.heatmap(
		x = nano.cnas,
		filename = BoutrosLab.utilities::generate.filename(fname.stem, 'cna_heatmap', 'png'),
		cluster.dimensions = clust.dim,
		clustering.method = 'diana',
		#clustering.method = 'centroid',
		xlab.label = 'Genes',
		ylab.cex = 2,
		ylab.label = 'Samples',
		xlab.cex = 2,
		rows.distance.method = dist.method,
		cols.distance.method = dist.method,
		colour.centering.value = 3,
		at = plot.at,
		colourkey.labels = if (rounded) { c(0,2,4) } else { plot.seq },
		colourkey.labels.at = if (rounded) { c(1,3,5) } else { plot.seq },
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = c('blue', 'white', 'red'),
		resolution = 600,
		...
		);
	}

