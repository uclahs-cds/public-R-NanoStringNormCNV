
# create a function for calling copy number changes
call.copy.number.state <- function (input, reference, per.chip = FALSE, chip.info = NULL, to.round = TRUE, to.log = FALSE, multi.factor = 2) {

	# assign header names
	header.names <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	# grep X and Y-chromosome genes
	x.genes <- input$Name[grep(x = input$Name, pattern = 'chrX')];
	y.genes <- input$Name[grep(x = input$Name, pattern = 'chrY')];

	# remove y.genes from input	  TO DO: get gender info to determine whether this should be done!
	#input <- input[!input$Name %in% y.genes,];

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

	# if user specified to.round then do the following
	if (to.round) {

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
			out.cna[,which.n] <- round(out.cna[,which.n, drop = FALSE], digits = 1);
			out.cna.round <- out.cna[, which.cna, drop = FALSE];
			}

		else { out.cna.round <- round(out.cna[,which.cna, drop = FALSE], digits = 1); }

		# loop over each gene
		for (row.ind in 1:nrow(out.cna.round)) {

			# get a list of genes with the following criteria (should 0.5 CN be classified as 1?)
			which.0 <- which(out.cna.round[row.ind,] <= 0.4);
			which.1 <- which(out.cna.round[row.ind,] >= 0.6 & out.cna.round[row.ind,] <= 1.4);
			which.2 <- which(out.cna.round[row.ind,] >= 1.6 & out.cna.round[row.ind,] <= 2.4);
			which.3 <- which(out.cna.round[row.ind,] >= 2.6 & out.cna.round[row.ind,] <= 3.4);
			which.4 <- which(out.cna.round[row.ind,] >= 3.6);

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
