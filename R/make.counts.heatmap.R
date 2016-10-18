make.counts.heatmap <- function(nano.counts, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL, clust.dim = 'both', clust.method = 'euclidean', print.ylab = NULL) {
	c.row <- covs.rows;
	c.col <- covs.cols;

	# set up plot labelling
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

	# check and set up covariates
	if (!is.null(c.row)) {
		# check if information is complete
		if (! all( colnames(nano.counts) %in% c.row$SampleID )) {
			stop("Must provide covariate information for every sample!");
			}
		# check if samples are ordered correctly
		if (! all( colnames(nano.counts) == c.row$SampleID )) {
			c.row <- c.row[match(colnames(nano.counts), c.row$SampleID),];
			}

		c.row <- c.row[, !(names(c.row) == 'SampleID'), drop = FALSE];
		rownames(c.row) <- NULL;

		# create row covariate object
		row.cov.obj <- NanoStringNormCNV::generate.plot.covariates(
			cov.info = c.row
			);
	} else {
		row.cov.obj <- NULL;
		}

	if (!is.null(c.col)) {
		# check if information is complete
		if (! all( rownames(nano.counts) %in% c.col$Name )) {
			stop("Must provide covariate information for every gene!");
			}
		# check if genes are ordered correctly
		if (! all( rownames(nano.counts) == c.col$Name )) {
			c.col <- c.col[match(rownames(nano.counts), c.col$Name),];
			rownames(c.col) <- NULL;
			}

		c.col <- c.col[, !(names(c.col) == 'Name'), drop = FALSE];
		rownames(c.col) <- NULL;

		# create column covariate object
		col.cov.obj <- NanoStringNormCNV::generate.plot.covariates(
			cov.info = c.col
			);
	} else {
		col.cov.obj <- NULL;
		}

	# set up legend
	if (!is.null(c.col) | !is.null(c.row)) {
		if (!is.null(c.col) & !is.null(c.row)) {
			cov.list <- mapply(c, list(c.col), list(c.row), SIMPLIFY = FALSE)[[1]];
		} else if (!is.null(c.col)) {
			cov.list <- as.list(c.col);
		} else if (!is.null(c.row)) {
			cov.list <- as.list(c.row);
			}
		covs.legend <- NanoStringNormCNV::generate.plot.legend(cov.info = cov.list);
	} else {
		covs.legend <- NULL;
		}

	# set up file name
	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# plot
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
		yaxis.tck = 0,
		yaxis.lab = print.ylab,
		yaxis.cex = 1,
		covariates = row.cov.obj,
		covariates.grid.border = list(col = 'black', lwd = 1.5),
		covariates.top = col.cov.obj,
		covariates.top.grid.border = list(col = 'black', lwd = 1.5),
		covariate.legends = covs.legend,
		legend.border.padding = 1.5,
		colour.scheme = c('white', 'black'),
		colourkey.cex = 2,
		colourkey.labels.at = key.labels.at,
		colourkey.labels = key.labels,
		at = seq(0, ceiling(max(nano.counts)), by = split.by),
		resolution = 600,
		axis.xlab.padding = 1.5,
		width = 7
		);
	}
