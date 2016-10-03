
make.cna.heatmap <- function(nano.cnas, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL, rounded = FALSE, clust.dim = 'both', centering.value = 2, ...) {

	# must add in random CNAs to make plot work if all values are identical
	if(length(unique(c(nano.cnas))) == 1){
		nano.cnas[1, 1] 							<- 1;
		nano.cnas[nrow(nano.cnas), ncol(nano.cnas)] <- 4;
		}

	# define parameters up front
	dist.method <- ifelse(rounded, 'jaccard', 'euclidean');
	plot.at 	<- seq(floor(min(nano.cnas, na.rm = TRUE)), ceiling(max(nano.cnas, na.rm = TRUE)), 0.1);
	plot.seq 	<- seq(min(plot.at), max(plot.at));

	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	BoutrosLab.plotting.general::create.heatmap(
		x = nano.cnas,
		filename = paste0(Sys.Date(), fname.stem, '_cna-heatmap.tiff'),
		cluster.dimensions = clust.dim,
		clustering.method = 'diana',
		xlab.label = 'Genes',
		ylab.cex = 2,
		ylab.label = 'Samples',
		xlab.cex = 2,
		rows.distance.method = dist.method,
		cols.distance.method = dist.method,
		at = plot.at,
		colourkey.labels = plot.seq,
		colourkey.labels.at = plot.seq,
		colour.centering.value = centering.value,
		colourkey.cex = 2,
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = c('blue', 'white', 'red'),
		resolution = 600,
		axis.xlab.padding = 1.5,
		...
		);
	}
