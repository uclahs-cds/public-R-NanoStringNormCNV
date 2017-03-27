make.invariant.probe.plot <- function(inv.probe.counts, tissue.type = NULL) {
	# process sample annotation
	if (!is.null(tissue.type)) {
		if (!all(c('SampleID', 'Type') %in% colnames(tissue.type))) {
			flog.warn("Missing tissue type information!");
			tissue.type <- NULL;
		} else {
			tissue.type <- tissue.type[, colnames(tissue.type) %in% c("SampleID", "Type")];
			}
		}

	# get mean counts per gene
	mean.counts <- apply(inv.probe.counts, 1, mean);

	# set up and create scatterplot
	splot.df <- melt(
		cbind(probe = c(paste0('probe', 1:nrow(inv.probe.counts))), inv.probe.counts),
		id.vars = 'probe'
		);
	
	colnames(splot.df) <- c("probe", "sample", "count");
	splot.df$cols <- 'black';
	splot.df$cols[splot.df$count < 100] <- 'red';
	splot.df$probe.id <- seq(1:nrow(inv.probe.counts));
	
	curves.list <- list();
	for (i in 1:nrow(inv.probe.counts)) {
		curve.expr <- function(x) {};
		body(curve.expr) <- bquote(log10(mean.counts[.(i)]));
		curves.list[[i]] <- curve.expr;
		}

	xaxis.range <- pretty(seq(1:nrow(inv.probe.counts)));
	xaxis.range <- xaxis.range[xaxis.range < nrow(inv.probe.counts)];

	BoutrosLab.plotting.general::create.scatterplot(
		formula = log10(count) ~ jitter(probe.id),
		data = splot.df[, c('count', 'probe.id')],
		filename = paste0(Sys.Date(), '_all-invariant-probe-counts_scatterplot.tiff'),
		main = 'Invariant Probe Counts',
		main.cex = 2,
		col = splot.df$cols,
		ylimits = c(1, max(log10(splot.df$count)) + 0.5),
		xlab.label = "Invariant Probe Number", 
		ylab.label = expression('log'[10]*"Count"),
		xlab.cex = 1.8,
		ylab.cex = 2,
		xaxis.lab = xaxis.range,
		xat = xaxis.range,
		xaxis.cex = 1.2,
		abline.h = 2,
		abline.col = 'red',
		alpha = 0.85,
		add.curves = TRUE,
		curves.from = seq(1:nrow(inv.probe.counts)) - 0.25,
		curves.to = seq(1:nrow(inv.probe.counts)) + 0.25,
		curves.col = 'cyan',
		curves.exprs = curves.list,
		width = 6 + 0.15 * nrow(inv.probe.counts),
		resolution = 500
		);

	# set up and create barplot, if applicable
	if (all(inv.probe.counts >= 100)) {
		flog.info("All invariant probe counts pass minimum threshold of 100");
	} else {
		# plot the number of counts < 100 per sample
		bplot.df <- data.frame(
			n.low.counts = apply(inv.probe.counts, 2, function(f) length(which(f < 100))),
			samples = colnames(inv.probe.counts)
			);
		bplot.df <- bplot.df[bplot.df$n.low.counts > 0 ,];
		bplot.df$sample.id <- 1:nrow(bplot.df);

		bar.cols   <- 'black';
		bar.legend <- NULL;
		if (!is.null(tissue.type)) {
			if (all(bplot.df$samples %in% tissue.type$SampleID)) {
				tissue.type <- tissue.type[tissue.type$SampleID %in% bplot.df$samples,];

				if (all(tolower(unique(as.character(tissue.type$Type))) %in% c('reference', 'tumour'))) {
					bar.cols <- as.character(tissue.type$Type);
					bar.cols[tolower(bar.cols) == 'tumour'] <- 'black';
					bar.cols[tolower(bar.cols) == 'reference'] <- 'red';

					bar.legend <- list(
						right = list(
							fun = draw.key,
							args = list(
								key = list(
									points = list(pch = 22, cex = 2, fill = c('black', 'red')),
									text   = list(lab = c('Tumour', 'Reference'))
									)
								)
							)
						);
					}
			} else {
				flog.warn("'Tissue type' sample IDs do not match 'invariant probe count' sample IDs!")
				}
			}
			
		BoutrosLab.plotting.general::create.barplot(
			formula = sample.id ~ n.low.counts,
			data = bplot.df[, c('sample.id', 'n.low.counts')],
			filename = paste0(Sys.Date(), '_low-invariant-probe-counts_barplot.tiff'),
			xlimits = c(0, nrow(inv.probe.counts) + 1),
			xlab.label = 'Invariant probes < 100 (N)',
			ylab.label = 'Sample IDs',
			xlab.cex = 1.6,
			yaxis.lab = bplot.df$samples,
			yat = seq(1, nrow(bplot.df), by = 1),
			yaxis.cex = 1,
			plot.horizontal = TRUE,
			legend = bar.legend,
			col = bar.cols,
			height = 5,
			width = 8,
			resolution = 600
			);
		}
	}
