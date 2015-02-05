
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

	out.cna <- get.tumour.normal.ratio(reference, n.chip, chip.info, out.cna);
	out.cna <- out.cna * multi.factor;

	# if user specified to round then do the following (based on NS recommendataions)
	if (thresh.method == 'round') {
		# segment using NS predifined thresholds
		out.cna.round <- apply.ns.cna.thresh(out.cna);

		# add the probe information back to out.cna.round
		out.cna.round <- cbind(
			input[,colnames(input)[colnames(input) %in% header.names], drop = FALSE],
			out.cna.round
			);

		# return out.cna.round
		return (out.cna.round);
		}
	else if(thresh.method == 'KD'){
		# segment using kernel density
		out.cna.round <- apply.kd.cna.thresh(out.cna, kd.vals);
		
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
