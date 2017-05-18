make.restriction.fragmentation.plot <- function(restr.data, low.ratio.samples = NULL){
	# reformat for plotting
	restr.df <- melt(data = t(restr.data));
	colnames(restr.df) <- unlist(strsplit("sample site count","\\s"));
	restr.df$sample.id <- rep(seq(1:ncol(restr.data)), 2);

	# barplot first
	bplot <- BoutrosLab.plotting.general::create.barplot(
		formula = count ~ sample.id,
		data = restr.df[, c('count', 'sample.id')],
		groups = restr.df$site,
		stack = TRUE,
		xlimits = c(0.5, col(restr.data) + 0.5),
		ylimits = c(0, max(restr.df$count) * 1.1),
		xaxis.tck = 0.5,
		border.col = default.colours(2),
		col = default.colours(2),
		add.rectangle = TRUE,
		col.rectangle = 'blue',
		alpha.rectangle = 0.5,
		xleft.rectangle = low.ratio.samples - 0.5,
		ybottom.rectangle = -1000,
		xright.rectangle = low.ratio.samples + 0.5,
		ytop.rectangle = max(restr.df$count) * 1.5,
		legend = list(
			  inside = list(
					fun = lattice::draw.key,
					args = list(
						key = list(
						   points = list(
								 col = 'black',
								 pch = 22,
								 cex = 1.6,
								 fill = default.colours(2)
								 ),
						   text = list(lab = c('A+B','C+D')),
						   padding.text = 2,
						   cex = 1
						   )
						),
				x = 0.95,
				y = 0.65
				)
			)
		);

	# scatter plot of ratio
	restr.data['ratio',] <- restr.data['C+D',] / restr.data['A+B',];
	splot.df <- data.frame(
		samples = colnames(restr.data),
		sample.id = seq(1:ncol(restr.data)),
		ratio = t(restr.data['ratio',]),
		cols = rep('black', ncol(restr.data)),
		stringsAsFactors = FALSE
		);
	
	if(! is.null(low.ratio.samples)){
		splot.df$cols[low.ratio.samples] <- 'red';
		}

	splot <- BoutrosLab.plotting.general::create.scatterplot(
		formula = log10(ratio) ~ sample.id,
		data = splot.df[, c('ratio', 'sample.id')],
		type = 'p',
		cex = 1,
		col = splot.df$cols,
		points.col = splot.df$cols,
		xlimits = c(0, col(restr.data)),
		abline.h = log10(10),
		abline.lty = 2
		);
	
	# combine plots into single figure
	pretty.bplot <- c(
		floor(head(pretty(restr.df$count), 1)),
		ceiling(tail(pretty(restr.df$count), 1))
		);
	pretty.splot <- c(
		floor(head(pretty(log10(splot.df$ratio)), 1)),
		ceiling(tail(pretty(log10(splot.df$ratio)), 1))
		);

	BoutrosLab.plotting.general::create.multiplot(
		plot.objects = list(bplot, splot),
		filename = paste0(Sys.Date(), '_restriction-fragmentation-ratios_multiplot.tiff'),
		main = 'Restriction Fragmentation Norm',
		main.cex = 2,
		panel.heights = c(1, 2.5),
		ylimits = list(
			c(pretty.bplot[1], pretty.bplot[2]),
			c(0, pretty.splot[2])
			# c(pretty.splot[1], pretty.splot[2])
			),
		xlab.label = 'Samples',
		ylab.label = c(expression('\t    log'[10]*'(ratio)'), '\t\tCount'),
		ylab.cex = 2,
		ylab.padding = 6,
		xaxis.cex = c(1, 0),
		yaxis.cex = c(1, 1),
		yat = list(
			pretty(restr.df$count),
			c(0:pretty.splot[2])
			# c(pretty.splot[1]:pretty.splot[2])
			),
		xaxis.alternating = 1,
		x.relation = 'same',
		y.relation = 'free',
		width = 20,
		resolution = 600
		);
	}
