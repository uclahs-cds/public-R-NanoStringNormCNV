apply.ns.cna.thresh <- function(ratio.data, cna.thresh = c(0.4, 1.5, 2.5, 3.5)) {
	# assign header names
	headers <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	# define sample columns
	which.cna <- colnames(ratio.data)[!colnames(ratio.data) %in% headers];
	which.n   <- which(colnames(ratio.data) %in% which.cna);

	# pull below into a separate object
	na.counts <- apply(
		X = ratio.data[, which.n, drop = FALSE],
		MARGIN = 2,
		FUN = function(f) { all(is.na(f)) }
		);

	# need to exclude columns with all NAs:
	# occurs when perchip = T and there are no reference samples on that chip!
	if (any(na.counts)) {
		all.na <- which(na.counts);
		print(paste("dropping:", all.na));
		which.n <- which.n[-all.na];
		which.cna <- which.cna[which.n];
		}
	cna.output <- ratio.data[, which.cna, drop = FALSE];

	# apply thresholds
	tmp.out <- ratio.data[ , which.cna, drop = FALSE];
	cna.output[tmp.out <= cna.thresh[1]] <- 0;
	cna.output[tmp.out > cna.thresh[1] & tmp.out <= cna.thresh[2]] <- 1;
	cna.output[tmp.out > cna.thresh[2] & tmp.out <= cna.thresh[3]] <- 2;
	cna.output[tmp.out > cna.thresh[3] & tmp.out <= cna.thresh[4]] <- 3;
	cna.output[tmp.out > cna.thresh[4]] <- 4;

	return(cna.output);
	}
