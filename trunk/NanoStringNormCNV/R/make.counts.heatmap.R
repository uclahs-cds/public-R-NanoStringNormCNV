make.counts.heatmap <- function(nano.counts, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL, clust.dim = 'both', clust.method = 'euclidean', print.ylab = NULL) {
	c.row <- covs.rows[covs.rows$SampleID %in% colnames(nano.counts),];
	c.col <- covs.cols[covs.cols$Name %in% rownames(nano.counts),];

	# set up plot labelling
	split.by  	  <- 1;
	key.labels.at <- NULL;
	key.labels 	  <- NULL;

	if (max(nano.counts, na.rm = TRUE) > 5000) {
		nano.counts <- log10(nano.counts + 1);
		split.by 	<- 0.1;
	} else if (max(nano.counts, na.rm = TRUE) == 1) {
		split.by 	  <- 0.5;
		key.labels.at <- c(0.25, 0.75);
		key.labels 	  <- c(0, 1);
		}

	# allow for single sample
	ylab.label <- 'Samples';
	if (ncol(nano.counts) == 1) {
		if (! is.null(print.ylab)) {
			print.ylab <- NULL;
			ylab.label <- paste0("Sample:\n", names(nano.counts));
			}
		nano.counts <- cbind(nano.counts, nano.counts);
		}

	# set up covariates and legend
	row.cov.obj <- NULL;
	col.cov.obj <- NULL;

	if (!is.null(c.col) | !is.null(c.row)) {
		# covariates
		cov.objs <- NanoStringNormCNV::generate.plot.covariates(
			plotting.data = nano.counts,
			sample.covariates = c.row,
			gene.covariates = c.col
			);
		row.cov.obj <- cov.objs[['sample']];
		col.cov.obj <- cov.objs[['gene']];

		# legend
		if (!is.null(c.col) & !is.null(c.row)) {
			cov.list <- mapply(c, list(c.col), list(c.row), SIMPLIFY = FALSE)[[1]];
		} else if (!is.null(c.col)) {
			cov.list <- as.list(c.col);
		} else if (!is.null(c.row)) {
			cov.list <- as.list(c.row);
			}
			
		cov.list <- cov.list[!(names(cov.list) %in% c('SampleID', 'Name'))];
		covs.legend <- NanoStringNormCNV::generate.plot.legend(cov.info = cov.list);
	} else {
		covs.legend <- NULL;
		}

	# set up file name
	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# plot
	BoutrosLab.plotting.general::create.heatmap(
		x = nano.counts,
		filename = paste0(Sys.Date(), fname.stem, '_counts-heatmap.tiff'),
		cluster.dimensions = clust.dim,
		rows.distance.method = clust.method,
		cols.distance.method = clust.method,
		xlab.label = 'Genes',
		ylab.label = ylab.label,
		xlab.cex = 2,
		ylab.cex = 2,
		yaxis.tck = 0,
		yaxis.lab = print.ylab,
		yaxis.cex = 1,
		covariates = row.cov.obj,
		covariates.grid.border = list(col = 'black', lwd = 1.5),
		covariates.top = col.cov.obj,
		covariates.top.grid.border = list(col = 'black', lwd = 1.5),
		legend.title.just = 'left',
		covariate.legends = covs.legend,
		legend.border.padding = 1.5,
		colour.scheme = c('white', 'black'),
		colourkey.cex = 2,
		colourkey.labels.at = key.labels.at,
		colourkey.labels = key.labels,
		at = seq(
			floor(min(nano.counts, na.rm = TRUE)),
			ceiling(max(nano.counts, na.rm = TRUE)),
			by = split.by
			),
		resolution = 600,
		axis.xlab.padding = 1.5,
		width = 7
		);
	}
