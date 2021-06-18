positive.control.qc <- function(raw.data) {
	# positive control target concentrations according to the NS nCounter CNV manual (pg 7)
	target.conc <- c(128, 32, 8, 2, 0.5, 0.125);
	target.conc <- target.conc * 171.23 + 214.12;

	# only keep the rows containing positive controls
	nano.pos <- raw.data['Positive' == raw.data$CodeClass,];
	nano.pos <- nano.pos[with(nano.pos, order(Name)),];

	# initiate object to store R squared values
	R2.results <- matrix(
		ncol = ncol(nano.pos) - 3,
		nrow = 1,
		dimnames = list("R2", colnames(nano.pos)[4:ncol(nano.pos)])
        );
	R2.results <- as.data.frame(R2.results);

	# for each sample, calculate R squared
	for (col.ind in 4:ncol(nano.pos)) {
		R2.results[1, colnames(nano.pos)[col.ind] == colnames(R2.results)] <- cor(
			nano.pos[,col.ind], target.conc
			)^2;
		}

	# throw a warning if all positive controls do not have R squared > 0.95
	if (!all(R2.results[1,] >= 0.95)) { warning('Not all probes have R^2 > 0.95'); }
	
	return(R2.results);
	}
