
make.counts.heatmap <- function(nano.counts, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL, clust.dim = 'both', clust.method = 'euclidean', print.ylab = NULL) {
	split.by  	  <- 1;
	key.labels.at <- NULL;
	key.labels 	  <- NULL;

	if (max(nano.counts) > 5000) {
		nano.counts <- log10(nano.counts + 1);
		split.by 	<- 0.1;
	} else if (max(nano.counts) == 1) {
		split.by 	  <- 0.5;
		key.labels.at <- c(0.25, 0.75);
		key.labels 	  <- c(0, 1);
		}

	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	BoutrosLab.plotting.general::create.heatmap(
		x = nano.counts,
		filename = paste0(Sys.Date(), fname.stem, '_counts_heatmap.tiff'),
		cluster.dimensions = clust.dim,
		rows.distance.method = clust.method,
		cols.distance.method = clust.method,
		xlab.label = 'Genes',
		ylab.label = 'Samples',
		xlab.cex = 2,
		ylab.cex = 2,
		yaxis.lab = print.ylab,
		yaxis.cex = 1,
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = c('white', 'black'),
		colourkey.cex = 2,
		colourkey.labels.at = key.labels.at,
		colourkey.labels = key.labels,
		at = seq(0, max(nano.counts), by = split.by),
		resolution = 600,
		axis.xlab.padding = 1.5,
		width = 7
		);
	}
