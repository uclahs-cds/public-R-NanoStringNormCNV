
# create a function for calling copy number changes
call.copy.number.state <- function (input, reference, per.chip = FALSE, chip.info = NULL, thresh.method = 'round', to.log = FALSE, multi.factor = 2, kd.vals = c(0.85,0.95)) {

	# Check input
	if(! thresh.method %in% qw("round KD kd none")){
		stop("Sorry method isn't currently supported. Please try one of round, KD, or none.");
		}

	if(toupper(thresh.method) == 'KD' & (2 != length(kd.vals) & 4 != length(kd.vals))){
		stop("Please specify two or four values for KD thresholds. The first should be for heterozygous and the second for homozygous if length 2. If length 4, the order should be hom deletion, het deletion, het gain, hom gain.");
		}

	# make sure kd values make sense
	if(toupper(thresh.method) == 'KD' & 2 == length(kd.vals) & kd.vals[1] > kd.vals[2]){
		stop("Invalid KD thresholds-- the first should be for heterozygous and the second for homozygous.");
		}
	if(toupper(thresh.method) == 'KD' & 4 == length(kd.vals) & (kd.vals[1] > kd.vals[2] | kd.vals[3] > kd.vals[4])){
		print(kd.vals);
		stop("Invalid KD thresholds-- the order should be hom deletion, het deletion, het gain, hom gain.");
		}

	# grep X and Y-chromosome genes
	x.genes <- input$Name[grep(x = input$Name, pattern = 'chrX')];
	y.genes <- input$Name[grep(x = input$Name, pattern = 'chrY')];
	if(length(x.genes) > 0){
		warning("*** AT THE MOMENT WE IGNORE GENES/PROBES ON CHRX AND CHRY!!!****");
		# remove x.genes from input	  TO DO: get gender info to determine whether this should be done!
		input <- input[!input$Name %in% x.genes,];
		}
	if(length(y.genes) > 0){
		warning("*** AT THE MOMENT WE IGNORE GENES/PROBES ON CHRX AND CHRY!!!****");
		# remove y.genes from input	  TO DO: get gender info to determine whether this should be done!
		input <- input[!input$Name %in% y.genes,];
		}

	### Analysis
	# assign header names
	header.names <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	# create empty data-frame to store data
	out.cna <- input;
	to.keep <- colnames(input);

	# if header is present remove headers from to.keep
	if (any(header.names %in% to.keep)) { to.keep <- to.keep[!to.keep %in% header.names]; }

	# make out.cna all NAs for data storage
	out.cna <- out.cna[,to.keep];
	out.cna[,to.keep] <- NA;

	# remove to.keep
	rm(to.keep);

	# see if user asks for per.chip
	if (per.chip)       { n.chip <- length(unique(chip.info$Chip)); }
	else if (!per.chip) { n.chip <- 1;                              }

	# define samples.to.loop here so the reference sample is placed at the end
	samples.to.loop <- c(colnames(out.cna)[!colnames(out.cna) %in% reference], reference);

	# loop over n.chip
	for (perm.i in 1:n.chip) {

		# define a tmp.ref every time before the start of a new iteration
		tmp.ref <- reference;

		# if per.chip is TRUE
		if (per.chip) {

			# get the tmp.ref for chip specifically
			this.chip <- paste0('Chip ', perm.i);
			tmp.ref <- chip.info$SampleID[chip.info$Chip %in% this.chip & chip.info$SampleID %in% reference];

			# when per chip is requested, but there are no ref samples on the chip, use pooled refs then
			if (length(tmp.ref) < 1) { next; }

			samples.to.loop <- chip.info$SampleID[this.chip == chip.info$Chip];
			samples.to.loop <- c(samples.to.loop[!samples.to.loop %in% tmp.ref], tmp.ref);
			}

		# if length of tmp.ref is greater than 1, take the average of the tmp.ref
		if (length(tmp.ref) > 1) {

			# take the average and create a new column in input as avg.ref
			input$avg.ref <- apply(X = input[,tmp.ref], MARGIN = 1, FUN = mean);

			# and change tmp.ref to 'avg.ref'
			tmp.ref <- 'avg.ref';
			}

		# loop over each sample first and then do the processing for reference samples!
		for (this.sample in samples.to.loop) {

			# divide each test samples probe value by corresponding probes in the ref samples
			if (any(0 == input[,tmp.ref])) { input[,tmp.ref][0 == input[,tmp.ref]] <- 1; }
			tmp.ratio <- input[,this.sample] / input[,tmp.ref];

			# multiply ratio by multiplication factor
			out.cna[,this.sample] <- tmp.ratio * multi.factor;
			}
		}

	# if user specified to round then do the following (based on NS recommendataions)
	if (thresh.method == 'round') {

		# round the numbers to the nearest integer ** rounding 1.5 and 2.5 to 2 and 3 respectively **
		which.cna <- colnames(out.cna)[!colnames(out.cna) %in% header.names];
		which.n   <- which(colnames(out.cna) %in% which.cna);

		# pull below into a separate object
		na.counts <- apply(
			X = out.cna[,which.n, drop = FALSE],
			MARGIN = 2,
			FUN = function(f) { all(is.na(f)) }
			);

		# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
		if (any(na.counts)) {
			all.na <- which(na.counts);
			which.n <- which.n[-all.na];
			print(paste("dropping:", all.na));
			which.cna <- which.cna[which.n];
			out.cna[,which.n] <- round(out.cna[,which.n, drop = FALSE], digits = 1);
			out.cna.round <- out.cna[, which.cna, drop = FALSE];
			}

		else { out.cna.round <- round(out.cna[,which.cna, drop = FALSE], digits = 1); }

		# loop over each gene
		for (row.ind in 1:nrow(out.cna.round)) {

			# get a list of genes with the following criteria (should 0.5 CN be classified as 1?)
			which.0 <- which(out.cna.round[row.ind,] <= 0.4);
			which.1 <- which(out.cna.round[row.ind,] > 0.4 & out.cna.round[row.ind,] <= 1.5);	# changed from 1.4 to 1.5
			which.2 <- which(out.cna.round[row.ind,] > 1.5 & out.cna.round[row.ind,] <= 2.5);	# changed from 2.4 to 2.5
			which.3 <- which(out.cna.round[row.ind,] > 2.5 & out.cna.round[row.ind,] <= 3.5);	# changed from 3.4 to 3.5
			which.4 <- which(out.cna.round[row.ind,] > 3.5);	# changed from >= 3.6 to >3.6

			# convert value to the nearest integer
			out.cna.round[row.ind, which.0] <- 0;
			out.cna.round[row.ind, which.1] <- 1;
			out.cna.round[row.ind, which.2] <- 2;
			out.cna.round[row.ind, which.3] <- 3;
			out.cna.round[row.ind, which.4] <- 4;
			}

		# add the probe information back to out.cna.round
		out.cna.round <- cbind(
			input[,colnames(input)[colnames(input) %in% header.names], drop = FALSE],
			out.cna.round
			);

		# return out.cna.round
		return (out.cna.round);
		}
	else if(thresh.method == 'KD'){
		which.cna <- colnames(out.cna)[!colnames(out.cna) %in% header.names];
		which.n   <- which(colnames(out.cna) %in% which.cna);

		# pull below into a separate object
		# pull below into a separate object
		na.counts <- apply(
			X = out.cna[,which.n, drop = FALSE],
			MARGIN = 2,
			FUN = function(f) { all(is.na(f)) }
			);

		# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
		if (any(na.counts)) {
			all.na <- which(na.counts);
			print(paste("dropping:", all.na));
			which.n <- which.n[-all.na];
			which.cna <- which.cna[which.n];
			out.cna.round <- out.cna[, which.cna, drop = FALSE];
			}
		else{
			out.cna.round <- out.cna[, which.cna, drop = FALSE];
			}

		# loop over each sample
		for (col.ind in 1:ncol(out.cna.round)) {
			if(2 == length(kd.vals)){
				cna.thresh.single <- get.sample.specific.cna.thresholds(method = 4, data = out.cna.round[ , col.ind], percent = kd.vals[1])[1:2];
				cna.thresh.multi <- get.sample.specific.cna.thresholds(method = 4, data = out.cna.round[ , col.ind], percent = kd.vals[2])[1:2];
				out.cna.round[ , col.ind] <- as.vector(call.cna.states(data.frame(log2ratio=out.cna.round[ , col.ind]), c(cna.thresh.multi[1], cna.thresh.single, cna.thresh.multi[2]))$CN) + 2;
			}else if(4 == length(kd.vals)){
				thresh <- vector(length = 4);
				thresh[1] <- get.sample.specific.cna.thresholds(method = 4, data = out.cna.round[ , col.ind], percent = kd.vals[1])[1];	# hom del
				thresh[2] <- get.sample.specific.cna.thresholds(method = 4, data = out.cna.round[ , col.ind], percent = kd.vals[2])[1];	# het del
				thresh[3] <- get.sample.specific.cna.thresholds(method = 4, data = out.cna.round[ , col.ind], percent = kd.vals[3])[2];	# het gain
				thresh[4] <- get.sample.specific.cna.thresholds(method = 4, data = out.cna.round[ , col.ind], percent = kd.vals[4])[2];	# hom gain
				out.cna.round[ , col.ind] <- as.vector(call.cna.states(data.frame(log2ratio=out.cna.round[ , col.ind]), thresh)$CN) + 2;
				}
			}
		# add the probe information back to out.cna.round
		out.cna.round <- cbind(
			input[,colnames(input)[colnames(input) %in% header.names], drop = FALSE],
			out.cna.round
			);
		return(out.cna.round);

		}

	# else return out.cna
	else {

		return(
			cbind(
				input[,colnames(input)[colnames(input) %in% header.names], drop = FALSE],
				out.cna
				)
			);
		}
	}
