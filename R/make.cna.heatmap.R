make.cna.heatmap <- function(nano.cnas, fname.stem = NULL, covs.rows = NULL, covs.cols = NULL, rounded = FALSE, clust.dim = 'both', centering.value = 2, min.cn = NULL, ...) {
	c.row <- covs.rows[covs.rows$SampleID %in% colnames(nano.cnas),];
	c.col <- covs.cols[covs.cols$Name %in% rownames(nano.cnas),];

	# must add in random CNAs to make plot work if all values are identical
	if (length(unique(c(nano.cnas))) == 1) {
		nano.cnas[1, 1] 							<- 0;
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

	if (nrow(nano.cnas) > 1) {
		# set up covariates and legend
		row.cov.obj <- NULL;
		col.cov.obj <- NULL;

		if (!is.null(c.col) | !is.null(c.row)) {
			# remove 'Type' if all samples are either 'Tumour' or 'Reference'
			if ("Type" %in% colnames(c.row)) {
				if (nlevels(as.factor(c.row$Type)) == 1) {
					c.row <- c.row[, !colnames(c.row) %in% "Type", drop = FALSE];
					}
				}

			# create covariates
			cov.objs <- NanoStringNormCNV::generate.plot.covariates(
				plotting.data = nano.cnas,
				sample.covariates = c.row,
				gene.covariates = c.col
				);
			row.cov.obj <- cov.objs[['sample']];
			col.cov.obj <- cov.objs[['gene']];

			# create legend
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

		# define parameters up front
		dist.method <- ifelse(rounded, 'jaccard', 'euclidean');
		if (is.null(min.cn)) { min.cn <- floor(min(nano.cnas, na.rm = TRUE)); }
		plot.at  <- seq(min.cn, ceiling(max(nano.cnas, na.rm = TRUE)), 0.1);
		plot.seq <- seq(min(plot.at), max(plot.at));

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
	} else {
		flog.warn("Unable to plot heatmap: requires more data!");
		}
	}
