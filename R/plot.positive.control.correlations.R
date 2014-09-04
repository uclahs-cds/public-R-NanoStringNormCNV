
plot.positive.control.correlations <- function(corr.df, fname.stem){		### TO DO: Add in covariates and allow clustering
	# order matrix according to correlation
	corr.df <- corr.df[order(corr.df[,1]) ,];
	create.heatmap(
		x = as.matrix(corr.df),
		filename = generate.filename(fname.stem, 'full_range', 'png'),
		main = 'Positive Probe Correlations',
		main.cex = 2,
		cluster.dimension = 'none',
		height = 2,
		scale.data=F,
		xat = seq(0, nrow(corr.df), by = 20),
		xaxis.lab = rownames(corr.df)[seq(20,nrow(corr.df), by = 20)],
		xlab.label = 'Samples',
		at = seq(0,1,by =0.01),
		colourkey.cex = 1.5,
		resolution = 600
		);
	create.heatmap(
		x = as.matrix(corr.df),
		filename = generate.filename(fname.stem, 'zoomed', 'png'),
		main = 'Positive Probe Correlations',
		main.cex = 2,
		cluster.dimension = 'none',
		height = 2,
		scale.data=F,
		xat = seq(0, nrow(corr.df), by = 20),
		xaxis.lab = rownames(corr.df)[seq(20,nrow(corr.df), by = 20)],
		xlab.label = 'Samples',
		colourkey.cex = 1.5,
		resolution = 600
		);
	}
