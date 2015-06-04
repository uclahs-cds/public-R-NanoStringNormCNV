

make.restr.digest.plot <- function(restr.data, low.ratio.samples = NULL){
	# reformat for plotting
	restr.df <- melt(data = t(restr.data));
	colnames(restr.df) <- unlist(strsplit("sample site count","\\s"));
	restr.df$sample.id <- rep(seq(1:ncol(restr.data)), 2);

	# barplot first
	bplot <- create.barplot(
		formula = count ~ sample.id,
		groups = restr.df$site,
		data = restr.df,
		stack = TRUE,
		ylab.label = 'Count',
		ylab.cex = 1.8,
		xaxis.cex = 1,
		xat = seq(20,nrow(restr.df), by=20),
		xaxis.tck = 0.5,
		border.col = default.colours(2),
		col = default.colours(2),
		ylimits = c(0, max(restr.df$count) * 1.1),
		xlimits = c(1,col(restr.data)),
		add.rectangle = TRUE,
		xleft.rectangle = low.ratio.samples - 0.5,
		ybottom.rectangle = -1000,
		xright.rectangle = low.ratio.samples + 0.5,
		ytop.rectangle = max(restr.df$count) * 1.5,
		col.rectangle = 'pink',
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
				x = 0.75,
				y = 0.50
				)
			),
	    width = 20,
		resolution = 600
		);

	# scatter plot of ratio
	restr.data['ratio', ] <- restr.data['C+D',]/restr.data['A+B',];
	splot.df <- data.frame(
		samples = colnames(restr.data),
		sample.id = seq(1:ncol(restr.data)),
		ratio = t(restr.data['ratio' , ]),
		cols = rep('black', ncol(restr.data))
		);
	
	if(! is.null(low.ratio.samples)){
		splot.df$cols[low.ratio.samples] <- 'red';
		}
	splot <- create.scatterplot(
		formula = log10(ratio) ~ sample.id,
		data = splot.df,
		type = 'p',
		points.col = splot.df$cols,
		col = splot.df$cols,
		cex = 1,
		xaxis.tck = 0.5,
		xlimits = c(1,col(restr.data)),
		abline.h = log10(50),	#1,
		ylab.label = 'C+D / A+B',
		ylab.cex = 2,
		xaxis.lab = seq(20, nrow(restr.df),by=20),
		xaxis.cex = 1,
		xat = seq(20,nrow(restr.df), by=20),
		width = 20,
		resolution = 600
		);
	
	create.multiplot(
		plot.objects = list(bplot,splot),
		main = 'Restriction Digestion Norm',
		filename = BoutrosLab.utilities::generate.filename('NanoString', 'restriction_digestion_ratios-mplot', 'png'),
		panel.heights = c(1, 3),
		main.cex = 2,
		xlab.label = 'Samples',
		x.relation = 'same',
		ylimits = list(c(0,5000), c(0,3)),
		xaxis.alternating = 1,
		merge.legends = TRUE,
		ylab.label = c(expression('\t\tlog'[10]*'(ratio)'), 'Count\t'),
		ylab.cex = 2,
		yaxis.cex = c(1, 1),
		retrieve.plot.labels = TRUE,
		width = 20,
		resolution = 600
		);
	}
