
apply.ns.cna.thresh <- function(tmr2ref){
	# assign header names
	headers <- c('Code.Class', 'CodeClass', 'Name', 'Accession');
	
	# round the numbers to the nearest integer ** rounding 1.5 and 2.5 to 2 and 3 respectively **
	which.cna <- colnames(tmr2ref)[!colnames(tmr2ref) %in% headers];
	which.n   <- which(colnames(tmr2ref) %in% which.cna);

	# pull below into a separate object
	na.counts <- apply(
		X = tmr2ref[,which.n, drop = FALSE],
		MARGIN = 2,
		FUN = function(f) { all(is.na(f)) }
		);

	# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
	if (any(na.counts)) {
		all.na <- which(na.counts);
		which.n <- which.n[-all.na];
		print(paste("dropping:", all.na));
		which.cna <- which.cna[which.n];
		tmr2ref[,which.n] <- round(tmr2ref[,which.n, drop = FALSE], digits = 1);
		cna.output <- tmr2ref[, which.cna, drop = FALSE];
		}

	else { cna.output <- round(tmr2ref[,which.cna, drop = FALSE], digits = 1); }

	# loop over each gene
	for (row.ind in 1:nrow(cna.output)) {

		# get a list of genes with the following criteria (should 0.5 CN be classified as 1?)
		which.0 <- which(cna.output[row.ind,] <= 0.4);
		which.1 <- which(cna.output[row.ind,] > 0.4 & cna.output[row.ind,] <= 1.5);	# changed from 1.4 to 1.5
		which.2 <- which(cna.output[row.ind,] > 1.5 & cna.output[row.ind,] <= 2.5);	# changed from 2.4 to 2.5
#		which.3 <- which(cna.output[row.ind,] > 2.7 & cna.output[row.ind,] <= 3.5);	# changed from 3.4 to 3.5	# emilie trying stricter threshold for gains
		which.3 <- which(cna.output[row.ind,] > 2.5 & cna.output[row.ind,] <= 3.5);	# changed from 3.4 to 3.5
		which.4 <- which(cna.output[row.ind,] > 3.5);	# changed from >= 3.6 to >3.6

		# convert value to the nearest integer
		cna.output[row.ind, which.0] <- 0;
		cna.output[row.ind, which.1] <- 1;
		cna.output[row.ind, which.2] <- 2;
		cna.output[row.ind, which.3] <- 3;
		cna.output[row.ind, which.4] <- 4;
		}

	return(cna.output);
	}
