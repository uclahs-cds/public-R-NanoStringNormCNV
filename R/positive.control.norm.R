
positive.control.norm <- function(nano.df){
	#--- positive control correlation ----------------------------------------------------------------#
	target.conc <- c(128, 32, 8, 2, 0.5, 0.125);	# According to the NS nCounter CNV manaual pg 7
	target.conc <- target.conc*171.23 + 214.12;

	# only keep the rows pertaining to positive code class
	nano.pos <- nano.df['Positive' == nano.df$CodeClass,];
	nano.pos <- nano.pos[with(nano.pos, order(Name)),];

	# initiate object to store R2
	R2.results <- data.frame(
        R2 = NA,
        SampleName = colnames(nano.pos)[4:ncol(nano.pos)],
        stringsAsFactors = FALSE
        );

	# for each sample, calculate R^2 and create scatterplot
	for (col.ind in 4:ncol(nano.pos)) {
		# get R^2 measure
		R2.results$R2[colnames(nano.pos)[col.ind] == R2.results$SampleName] <- cor(
			nano.pos[,col.ind], target.conc
			)^2;
		}

	# throw a warning if not all positive controls have R^2 > 0.95
	if (!all(R2.results$R2 >= 0.95)) { warning('Not all probes have R^2 > 0.95'); }
	
	return(R2.results);
	}
