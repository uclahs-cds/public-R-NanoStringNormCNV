
# plot the counts for the invariant probes. It is BAD to have counts < 100 for these probes as it signifies
# low DNA input, leading to unreliable CNA calls. Especially if the low counts are in the reference samples!!

make.invariant.probe.plot <- function(inv.probe.counts, fname.stem){
	mean.counts <- apply(inv.probe.counts, 1, mean);
	splot.df <- melt(cbind(probe=c(paste0('probe', 1:nrow(inv.probe.counts))), inv.probe.counts));
	colnames(splot.df) <- qw("probe sample count");
	splot.df$cols <- 'black';
	splot.df$cols[splot.df$count < 100] <- 'red';
	splot.df$probe.id <- seq(1:nrow(inv.probe.counts));
	
	curves.list <- list(
		function(x) log10(mean.counts[1]),
		function(x) log10(mean.counts[2]),
		function(x) log10(mean.counts[3]),
		function(x) log10(mean.counts[4]),
		function(x) log10(mean.counts[5]),
		function(x) log10(mean.counts[6]),
		function(x) log10(mean.counts[7]),
		function(x) log10(mean.counts[8]),
		function(x) log10(mean.counts[9])
#		function(x) log10(mean.counts[10])
		);

	create.scatterplot(
		log10(count) ~ jitter(probe.id),
		data = splot.df,
		filename = generate.filename(fname.stem, 'counts', 'png'),
		main = 'Invariant Probe Counts',
		main.cex = 2,
		col = splot.df$cols,
		xlab.label = "Invariant Probe Number", 
		xat = seq(1:nrow(inv.probe.counts)),
		xlab.cex = 1.8,
		xaxis.cex = 1.2,
		xaxis.lab = seq(1:nrow(inv.probe.counts)),
		ylab.label = expression('log'[10]*"Count"),
		abline.h = 2,
		abline.col = 'red',
		ylab.cex = 2,
		ylimits = c(1, max(log10(splot.df$count))+0.5),
		alpha = 0.85,
		add.curves = TRUE,
		curves.from = (seq(1:nrow(inv.probe.counts)) - 0.45),
		curves.to = (seq(1:nrow(inv.probe.counts)) + 0.45),
		curves.col = 'cyan',
		curves.exprs = curves.list,
		resolution = 500
		);


	# TO DO: colour bars accoring to sample type
	# plot the # of counts < 100 per sample
	bplot.df <- data.frame(
		n.low.counts = apply(inv.probe.counts, 2, function(f) length(which(f<100))),
		samples = colnames(inv.probe.counts)
		);
	bplot.df <- bplot.df[bplot.df$n.low.counts > 0 ,];
	bplot.df$sample.id <- seq(1:nrow(bplot.df));
	create.barplot(
		sample.id ~ n.low.counts,# ~ sample.id,
		data = bplot.df,
		filename = generate.filename(fname.stem, 'low_counts_per_sample', 'png'),
		yaxis.lab = bplot.df$samples,
		yat = seq(1,nrow(bplot.df), by = 1),
		yaxis.cex = 1,
		xlab.label = 'Invariant probes < 100 (N)',
		xlab.cex = 1.6,
		xlimits = c(0, nrow(inv.probe.counts) + 1),
		xat = seq(0, nrow(inv.probe.counts), by= 2),
		plot.horizontal = TRUE,
		resolution = 600
		);
	}
