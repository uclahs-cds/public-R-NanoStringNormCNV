make.positive.control.plot <- function(correlations, covs = NULL, print.x.labels = TRUE) {
	# order values by size
	correlations <- correlations[, order(correlations[1,])];

	# set up covariates and legend
	if (! 'SampleID' %in% colnames(covs)) {
		flog.warn("Unable to plot covariates with sample IDs!");
		covs <- NULL;
	} else {
		covs <- covs[, colnames(covs) %in% c('SampleID', 'Type', 'Cartridge')];
		}

	if (!is.null(covs) & ncol(covs) > 1) {
		# covariates
		cov.objs <- NanoStringNormCNV:::generate.plot.covariates(
			plotting.data = correlations,
			sample.covariates = covs
			);
		cov.obj <- cov.objs[['sample']];
		
		# legend
		covs <- covs[, names(covs) != 'SampleID'];
		covs.legend <- NanoStringNormCNV:::generate.plot.legend(cov.info = as.list(covs));

		clust.dim <- 'columns';
	} else {
		covs.legend <- NULL;
		cov.obj <- NULL;
		clust.dim <- 'none';
		}

	# set up plot parameters
	if (print.x.labels) {
		xlab.label <- 'Sample ID';
		xaxis.labels <- colnames(correlations);
		xat <- 1:length(xaxis.labels);
		plot.height <- 4;
	} else {
		xaxis.labels <- xat <- NULL;
		xlab.label <- '';
		plot.height <- 3;
		}

	# plot
	BoutrosLab.plotting.general::create.heatmap(
		x = t(rbind(correlations, correlations)),
		filename = paste0(Sys.Date(), '_positive-control-correlations_full-range.tiff'),
		cluster.dimensions = clust.dim,
		cols.distance.method = 'euclidean',
		main = 'Positive Probe Correlations',
		main.cex = 1.8,
		scale.data = FALSE,
		xlab.label = xlab.label,
		xlab.cex = 1.2,
		xaxis.lab = xaxis.labels,
		xat = xat,
		xaxis.cex = 0.5,
		yaxis.tck = 0,
		covariates.top.grid.border = list(col = 'black', lwd = 1.5),
		covariates.top = cov.obj,
		covariate.legends = covs.legend,
		legend.title.just = 'left',
		colour.scheme = c('red', 'white', 'blue'),
		colour.centering.value = 0.5,
		colourkey.cex = 1.5,
		at = seq(0, 1, by = 0.01),
		right.padding = 2,
		height = plot.height,
		width = 5 + ncol(correlations) * 0.09,
		resolution = 600
		);
	
	BoutrosLab.plotting.general::create.heatmap(
		x = t(rbind(correlations, correlations)),
		filename = paste0(Sys.Date(), '_positive-control-correlations_zoomed-in.tiff'),
		cluster.dimensions = clust.dim,
		cols.distance.method = 'euclidean',
		main = 'Positive Probe Correlations',
		main.cex = 1.8,
		scale.data = FALSE,
		xlab.label = xlab.label,
		xlab.cex = 1.2,
		xaxis.lab = xaxis.labels,
		xat = xat,
		xaxis.cex = 0.5,
		yaxis.tck = 0,
		covariates.top.grid.border = list(col = 'black', lwd = 1.5),
		covariates.top = cov.obj,
		covariate.legends = covs.legend,
		legend.title.just = 'left',
		colour.scheme = c('white','blue'),
		colourkey.cex = 1.5,
		right.padding = 2,
		height = plot.height,
		width = 5 + ncol(correlations) * 0.09,
		resolution = 600
		);
	}

