make.cna.heatmap <- function(nano.cnas, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL, rounded = FALSE, clust.dim = 'both', centering.value = 2, ...) {
	c.row <- covs.rows;
	c.col <- covs.cols;

	# must add in random CNAs to make plot work if all values are identical
	if(length(unique(c(nano.cnas))) == 1) {
		nano.cnas[1, 1] 							<- 1;
		nano.cnas[nrow(nano.cnas), ncol(nano.cnas)] <- 4;
		}

	# must remove probes that contain any NAs
	na.probes <- as.vector(which(is.na(rowSums(nano.cnas))));
	if (length(na.probes) > 0) {
		flog.warn(paste0(
			"Removing the following genes from plot due to NA values in one or more samples:\n",
			paste("\t", rownames(nano.cnas)[na.probes], collapse = "\n")
			));
		nano.cnas <- nano.cnas[-na.probes,];
		}

	# check and set up covariates
	if (!is.null(c.row)) {
		# check completeness and order
		if (! all( colnames(nano.cnas) %in% c.row$SampleID )) {
			stop("Must provide covariate information for every sample!");
		} else {
			c.row <- c.row[match(colnames(nano.cnas), c.row$SampleID),];
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
		# check completeness and order
		if (! all( rownames(nano.cnas) %in% c.col$Name )) {
			stop("Must provide covariate information for every gene!");
		} else {
			c.col <- c.col[match(rownames(nano.cnas), c.col$Name),];
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

	# define parameters up front
	dist.method <- ifelse(rounded, 'jaccard', 'euclidean');
	plot.at 	<- seq(floor(min(nano.cnas, na.rm = TRUE)), ceiling(max(nano.cnas, na.rm = TRUE)), 0.1);
	plot.seq 	<- seq(min(plot.at), max(plot.at));

	# set up file name
	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# plot
	BoutrosLab.plotting.general::create.heatmap(
		x = nano.cnas,
		filename = paste0(Sys.Date(), fname.stem, '_cna-heatmap.tiff'),
		cluster.dimensions = clust.dim,
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
		covariates.grid.border = list(col = 'black', lwd = 1.5),
		covariates = row.cov.obj,
		covariates.top.grid.border = list(col = 'black', lwd = 1.5),
		covariates.top = col.cov.obj,
		covariate.legends = covs.legend,
		colour.scheme = c('blue', 'white', 'red'),
		resolution = 600,
		axis.xlab.padding = 1.5,
		...
		);
	}
