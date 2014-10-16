
# create a function for calling copy number changes
call.copy.number.state <- function (input, reference, per.chip = FALSE, chip.info = NULL, to.round = TRUE, to.log = FALSE, multi.factor = 2) {

	# assign header names
	header.names <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	# grep X and Y-chromosome genes
	x.genes <- input$Name[grep(x = input$Name, pattern = 'chrX')];
	y.genes <- input$Name[grep(x = input$Name, pattern = 'chrY')];

	# remove y.genes from input	  TO DO: get gender info to determine whether this should be done!
	input <- input[!input$Name %in% y.genes,];

	# assess input file and create an empty data-frame for data storage
	if (any(header.names %in% colnames(input))) {

			# check if these three columns are in the header
			out.cna   <- input;
			which.cna <- colnames(input)[!colnames(input) %in% header.names];

			# make the rest of the data NA and only keeping the headers
			out.cna[,which.cna] <- NA;
			}
	else{
		out.cna <- input;
		out.cna <- NA;
		}

	# see if user asks for per.chip
	if (per.chip) {
			n.chip     <- length(unique(chip.info$Chip));
			ref.backup <- reference;
			}

	else if (!per.chip) { n.chip <- 1; }

	# loop over n.chip
	for (perm.i in 1:n.chip) {

			# if per.chip is TRUE
			if (per.chip) {

					# get the reference for chip specifically
					this.chip <- paste0('Chip ', perm.i);
					to.keep   <- chip.info$Chip %in% this.chip & chip.info$SampleID %in% ref.backup;
					reference <- chip.info$SampleID[to.keep];

					# when per chip is requested, but there are no reference samples on the chip, use pooled references then
					if(length(reference) < 1){
						next;
#						to.keep   <- chip.info$SampleID %in% ref.backup;
#						reference <- chip.info$SampleID[to.keep];
						}
					which.cna <- chip.info$SampleID[this.chip == chip.info$Chip];
					which.cna <- c(which.cna[!which.cna %in% reference], reference);
					}

			# if length of reference is greater than 1, take the average of the reference
			if (length(reference) > 1) {

					# take the average and create a new column in input as avg.ref
					input$avg.ref <- apply(X = input[,reference], MARGIN = 1, FUN = mean);

					# and change reference to 'avg.ref'
					reference <- 'avg.ref';
					}

			# loop over each sample
			for (this.sample in which.cna) {
					# divide each test samples probe value by corresponding probes in the reference samples
				if(any(input[,reference] == 0)){
						input[,reference][input[,reference] == 0] <- 1;
						}
					tmp.ratio <- input[,this.sample] / input[,reference];

					# multiply ratio by multiplication factor
					out.cna[,this.sample] <- tmp.ratio * multi.factor;
					}
			}

	# if user specified to.round then do the following
	if (to.round) {

			# round the numbers to the nearest integer ** rounding 1.5 and 2.5 to 2 and 3 respectively **
			which.cna     <- colnames(out.cna)[!colnames(out.cna) %in% header.names];
			which.n <- which(colnames(out.cna) %in% which.cna);
			# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
			if(any(apply(out.cna[, which.n, drop = FALSE], 2, function(f) all(is.na(f))))){
				all.na <- which(apply(out.cna[, which.n, drop = FALSE], 2, function(f) all(is.na(f))));
				which.n <- which.n[-all.na];
				out.cna[,which.n] <- round(out.cna[,which.n, drop = FALSE], digits = 1);
				out.cna.round <- out.cna[, which.cna, drop = FALSE];
			}else{
				out.cna.round <- round(out.cna[,which.cna, drop = FALSE], digits = 1);
				}

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
					out.cna[,colnames(input) %in% header.names, drop = FALSE],
					out.cna.round
					);

			# return out.cna.round
			return(out.cna.round);
			}

	# else return out.cna
	else { return(out.cna); }
	}
