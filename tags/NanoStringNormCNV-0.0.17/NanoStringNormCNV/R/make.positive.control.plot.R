
make.positive.control.plot <- function(corr.df, fname.stem){		### TO DO: Add in covariates and allow clustering
	# order matrix according to correlation
	corr.df <- corr.df[order(corr.df[,1]), 1, drop = FALSE];
	rownames(corr.df) <- NULL;
	
	BoutrosLab.plotting.general::create.heatmap(
		x = as.matrix(corr.df),
		filename = BoutrosLab.utilities::generate.filename(fname.stem, 'full_range', 'png'),
		main = 'Positive Probe Correlations',
		main.cex = 2,
		cluster.dimensions = 'none',
		height = 2,
		scale.data = FALSE,
		xat = seq(0, nrow(corr.df), by = 20),
		xaxis.lab = c(0, rownames(corr.df)[seq(20, nrow(corr.df), by = 20)]),
		xlab.label = 'Samples',
		at = seq(0, 1, by = 0.01),
		colourkey.cex = 1.5,
		resolution = 600,
		right.padding = 2
		);
	
	BoutrosLab.plotting.general::create.heatmap(
		x = as.matrix(corr.df),
		filename = BoutrosLab.utilities::generate.filename(fname.stem, 'zoomed', 'png'),
		main = 'Positive Probe Correlations',
		main.cex = 2,
		cluster.dimensions = 'none',
		height = 2,
		scale.data = FALSE,
		xat = seq(0, nrow(corr.df), by = 20),
		xaxis.lab = c(0, rownames(corr.df)[seq(20, nrow(corr.df), by = 20)]),
		xlab.label = 'Samples',
		colourkey.cex = 1.5,
		resolution = 600,
		right.padding = 2
		);
	}

