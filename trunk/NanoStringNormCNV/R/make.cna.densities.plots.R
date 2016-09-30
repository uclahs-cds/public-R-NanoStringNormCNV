
make.cna.densities.plots <- function(cnas, xlab, fname.stem = NULL) {

	# round values to 5 decimal places
	cnas <- round(cnas, digits = 5);

	if (!is.null(fname.stem)) { fname.stem <- paste0("_", fname.stem); }

	# per gene plots
	gene.list 		 <- lapply(seq(1:nrow(cnas)), function(f) cnas[f,]);
	names(gene.list) <- paste0('gene', 1:length(gene.list));

	BoutrosLab.plotting.general::create.densityplot(
		x = gene.list,
		filename = paste0(Sys.Date(), fname.stem, '_gene-densityplot.tiff'),
		xlab.label = xlab,
		xlab.cex = 2,
		main = 'Density per gene',
		main.cex = 2,
		add.rectangle = TRUE,
		xleft.rectangle = c(-5, 0.5, 1.5, 2.5, 3.5),
		ybottom.rectangle = c(0, 0, 0, 0, 0),
		xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 15),
		ytop.rectangle = c(50, 50, 50, 50, 50),
		col.rectangle = c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick'),
		alpha.rectangle = 0.6,
		lwd = 0.5,
		lty = 'dotted',
		resolution = 600
		);

	# per sample plots
	sample.list 	   <- lapply(seq(1:ncol(cnas)), function(f) cnas[, f]);
	names(sample.list) <- paste0('sample', 1:length(sample.list));

	BoutrosLab.plotting.general::create.densityplot(
		x = sample.list,
		filename = paste0(Sys.Date(), fname.stem, '_sample-densityplot.tiff'),
		xlab.label = xlab,
		xlab.cex = 2,
		main = 'Density per sample',
		main.cex = 2,
		add.rectangle = TRUE,
		xleft.rectangle = c(-5, 0.5, 1.5, 2.5, 3.5),
		ybottom.rectangle = c(0, 0, 0, 0, 0),
		xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 15),
		ytop.rectangle = c(50, 50, 50, 50, 50),
		col.rectangle = c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick'),
		alpha.rectangle = 0.6,
		lwd = 0.5,
		lty = 'dotted',
		resolution = 600
		);

	# per sample logged plots
	if (min(cnas) > 0) {
		sample.list 	   <- lapply(seq(1:ncol(cnas)), function(f) log2(cnas[, f]));
		names(sample.list) <- paste0('sample', 1:length(sample.list));

		BoutrosLab.plotting.general::create.densityplot(
			x = sample.list,
			filename = paste0(Sys.Date(), fname.stem, '_sample-logged-densityplot.tiff'),
			xlab.label = xlab,
			xlab.cex = 2,
			main = 'Density per sample (logged)',
			main.cex = 2,
			add.rectangle = TRUE,
			xleft.rectangle = log2(c(0.001, 0.5, 1.5, 2.5, 3.5)),
			ybottom.rectangle = c(0, 0, 0, 0, 0),
			xright.rectangle = log2(c(0.5, 1.5, 2.5, 3.5, 15)),
			ytop.rectangle = c(50, 50, 50, 50, 50),
			col.rectangle = c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick'),
			alpha.rectangle = 0.6,
			lwd = 0.5,
			lty = 'dotted',
			resolution = 600
			);
		}
	
	}
