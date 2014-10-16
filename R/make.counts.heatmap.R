
make.counts.heatmap <- function(nano.counts, fname.stem, covs.rows = NULL, covs.cols = NULL, covs.legend = NULL, clust.dim = 'both', clust.method='euclidean', print.ylab = NULL){
	split.by <-  1;
	if(max(nano.counts) > 5000){
		nano.counts <- log10(nano.counts + 1);
		split.by <- 0.1;
		}
	else if(max(nano.counts) ==1){
		split.by <- 0.5;
		}
	create.heatmap(
		x = nano.counts,
		filename = generate.filename(fname.stem, 'counts_heatmap', 'png'),
		cluster.dimensions = clust.dim,
		rows.distance.method = clust.method,
		cols.distance.method = clust.method,
		xlab.label = 'Genes',
		ylab.cex = 2,
		ylab.label = 'Samples',
		xlab.cex = 2,
		yaxis.lab = print.ylab,
		yaxis.cex = 1,
		covariates = covs.rows,
		covariates.top = covs.cols,
		covariate.legends = covs.legend,
		colour.scheme = c('white', 'black'),
		at = seq(0, max(nano.counts)+0.5, by = split.by),
		resolution = 600
		);

	}
