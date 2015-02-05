
apply.ns.cna.thresh <- function(output, headers){
	# round the numbers to the nearest integer ** rounding 1.5 and 2.5 to 2 and 3 respectively **
	which.cna <- colnames(output)[!colnames(output) %in% headers];
	which.n   <- which(colnames(output) %in% which.cna);

	# pull below into a separate object
	na.counts <- apply(
		X = output[,which.n, drop = FALSE],
		MARGIN = 2,
		FUN = function(f) { all(is.na(f)) }
		);

	# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
	if (any(na.counts)) {
		all.na <- which(na.counts);
		which.n <- which.n[-all.na];
		print(paste("dropping:", all.na));
		which.cna <- which.cna[which.n];
		output[,which.n] <- round(output[,which.n, drop = FALSE], digits = 1);
		output.round <- output[, which.cna, drop = FALSE];
		}

	else { output.round <- round(output[,which.cna, drop = FALSE], digits = 1); }

	# loop over each gene
	for (row.ind in 1:nrow(output.round)) {

		# get a list of genes with the following criteria (should 0.5 CN be classified as 1?)
		which.0 <- which(output.round[row.ind,] <= 0.4);
		which.1 <- which(output.round[row.ind,] > 0.4 & output.round[row.ind,] <= 1.5);	# changed from 1.4 to 1.5
		which.2 <- which(output.round[row.ind,] > 1.5 & output.round[row.ind,] <= 2.5);	# changed from 2.4 to 2.5
#		which.3 <- which(output.round[row.ind,] > 2.7 & output.round[row.ind,] <= 3.5);	# changed from 3.4 to 3.5	# emilie trying stricter threshold for gains
		which.3 <- which(output.round[row.ind,] > 2.5 & output.round[row.ind,] <= 3.5);	# changed from 3.4 to 3.5
		which.4 <- which(output.round[row.ind,] > 3.5);	# changed from >= 3.6 to >3.6

		# convert value to the nearest integer
		output.round[row.ind, which.0] <- 0;
		output.round[row.ind, which.1] <- 1;
		output.round[row.ind, which.2] <- 2;
		output.round[row.ind, which.3] <- 3;
		output.round[row.ind, which.4] <- 4;
		}

	return(output.round);
	}
