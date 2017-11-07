make.cna.densities.plots <- function(nano.cnas, fname.stem = NULL) {

	# check for minimum of 2 samples
	if (ncol(nano.cnas) < 2) {
		flog.error("Must have at least 2 samples to plot CNA densities!");
		}

	# round values to 5 decimal places
	nano.cnas <- round(nano.cnas, digits = 5);

	# must remove probes that contain any NAs
	na.probes <- as.vector(which(is.na(rowSums(nano.cnas))));
	if (length(na.probes) > 0) {
		flog.warn(paste0(
			"Removing the following genes from plot due to NA values in one or more samples:\n",
			paste("\t", rownames(nano.cnas)[na.probes], collapse = "\n")
			));
		nano.cnas <- nano.cnas[-na.probes, , drop = FALSE];
		}

	# set up file name
	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# set up plot colours
	plot.colours <- c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick');

	if (nrow(nano.cnas) > 0) {
		# plot per gene
		gene.list 		 <- lapply(seq(1:nrow(nano.cnas)), function(f) nano.cnas[f,]);
		names(gene.list) <- paste0('gene', 1:length(gene.list));

		BoutrosLab.plotting.general::create.densityplot(
			x = gene.list,
			filename = paste0(Sys.Date(), fname.stem, '_gene-densityplot.tiff'),
			xlab.label = expression('CNA'),
			xlab.cex = 2,
			main = 'Density per gene',
			main.cex = 2,
			add.rectangle = TRUE,
			xleft.rectangle = c(-100, 0.5, 1.5, 2.5, 3.5),
			ybottom.rectangle = rep(0, 5),
			xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 100),
			ytop.rectangle = rep(length(gene.list) + 3, 5),
			col.rectangle = plot.colours,
			alpha.rectangle = 0.6,
			lwd = 0.5,
			lty = 'dotted',
			resolution = 600
			);
	} else {
		flog.warn("Unable to plot gene density data: requires more data!");		
		}

	if (nrow(nano.cnas) > 1) {
		# plot per sample
		sample.list 	   <- lapply(seq(1:ncol(nano.cnas)), function(f) nano.cnas[, f]);
		names(sample.list) <- paste0('sample', 1:length(sample.list));

		BoutrosLab.plotting.general::create.densityplot(
			x = sample.list,
			filename = paste0(Sys.Date(), fname.stem, '_sample-densityplot.tiff'),
			xlab.label = expression('CNA'),
			xlab.cex = 2,
			main = 'Density per sample',
			main.cex = 2,
			add.rectangle = TRUE,
			xleft.rectangle = c(-100, 0.5, 1.5, 2.5, 3.5),
			ybottom.rectangle = rep(0, 5),
			xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 100),
			ytop.rectangle = rep(length(sample.list) + 3, 5),
			col.rectangle = plot.colours,
			alpha.rectangle = 0.6,
			lwd = 0.5,
			lty = 'dotted',
			resolution = 600
			);

		# plot per sample logged
		if (min(nano.cnas, na.rm = TRUE) > 0) {
			sample.list 	   <- lapply(seq(1:ncol(nano.cnas)), function(f) log2(nano.cnas[, f]));
			names(sample.list) <- paste0('sample', 1:length(sample.list));

			BoutrosLab.plotting.general::create.densityplot(
				x = sample.list,
				filename = paste0(Sys.Date(), fname.stem, '_sample-logged-densityplot.tiff'),
				xlab.label = expression('CNA'),
				xlab.cex = 2,
				main = 'Density per sample (logged)',
				main.cex = 2,
				add.rectangle = TRUE,
				xleft.rectangle = log2(c(0.0001, 0.5, 1.5, 2.5, 3.5)),
				ybottom.rectangle = rep(0, 5),
				xright.rectangle = log2(c(0.5, 1.5, 2.5, 3.5, 1000)),
				ytop.rectangle = rep(length(sample.list) + 3, 5),
				col.rectangle = plot.colours,
				alpha.rectangle = 0.6,
				lwd = 0.5,
				lty = 'dotted',
				resolution = 600
				);
			}
	} else {
		flog.warn("Unable to plot sample density data: requires more data!");
		}
	
	}
