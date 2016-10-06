
make.positive.control.plot <- function(correlations) {		### TO DO: Add in covariates and allow clustering
	# order correlation by size
	correlations <- sort(correlations);

	# set up x-axis
	xaxis.range <- pretty(seq(0, length(correlations)), n = 4);
	xaxis.range <- xaxis.range[xaxis.range < length(correlations)];
	
	BoutrosLab.plotting.general::create.heatmap(
		x = as.matrix(correlations),
		filename = paste0(Sys.Date(), '_positive-control-correlations_full-range.tiff'),
		cluster.dimensions = 'none',
		main = 'Positive Probe Correlations',
		main.cex = 1.8,
		scale.data = FALSE,
		xlab.label = 'Samples',
		xlab.cex = 1.2,
		xat = xaxis.range,
		xaxis.lab = xaxis.range,
		yaxis.tck = 0,
		colour.scheme = c('red', 'white', 'blue'),
		colour.centering.value = 0.5,
		colourkey.cex = 1.5,
		at = seq(0, 1, by = 0.01),
		right.padding = 2,
		height = 2,
		resolution = 600
		);
	
	BoutrosLab.plotting.general::create.heatmap(
		x = as.matrix(correlations),
		filename = paste0(Sys.Date(), '_positive-control-correlations_zoomed-in.tiff'),
		cluster.dimensions = 'none',
		main = 'Positive Probe Correlations',
		main.cex = 1.8,
		scale.data = FALSE,
		xlab.label = 'Samples',
		xlab.cex = 1.2,
		xat = xaxis.range,
		xaxis.lab = xaxis.range,
		yaxis.tck = 0,
		colour.scheme = c('white','blue'),
		colourkey.cex = 1.5,
		right.padding = 2,
		height = 2,
		resolution = 600
		);
	}

