
get.tumour.normal.ratio <- function(ns.counts, ref, chips.info, per.chip = FALSE){
	# create empty data-frame to store data
	output <- ns.counts;

	cols.to.keep   <- colnames(ns.counts);
	cols.to.remove <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	if (any(cols.to.remove %in% cols.to.keep)) {
		cols.to.keep <- cols.to.keep[!cols.to.keep %in% cols.to.remove];
		}
	if (! any(chips.info$SampleID %in% colnames(ns.counts))) {
		chips.info <- chips.info[match(colnames(ns.counts), chips.info$SampleID),];
		}

	output <- output[, cols.to.keep, drop = FALSE];
	output[, cols.to.keep] <- NA;

	# see if user asks for per.chip
	if (per.chip) {
		chips <- unique(chips.info$cartridge);
		if (length(chips) < 1) {
			flog.warn("Cannot process data per chip: missing cartridge (chip) information!");
			per.chip <- 0;
			}
		}
	if (!per.chip) {
		chips <- 'combined';
		}

	# define sample order here so the reference sample is placed at the end
	samples.to.loop <- c(colnames(output)[!colnames(output) %in% ref], ref);

	# loop over chips
	for (this.chip in chips) {

		# define a tmp.ref every time before the start of a new iteration
		tmp.ref <- ref;

		if (this.chip != 'combined') {

			# get the tmp.ref for given chip
			tmp.ref <- chips.info$SampleID[chips.info$cartridge %in% this.chip & chips.info$SampleID %in% ref];

			# skip when per chip is requested but there are no ref samples on the chip
			if (length(tmp.ref) < 1) { 
				flog.warn(paste0("No reference samples on cartridge ", this.chip, ". Try calling CNAs with per.chip = FALSE"));
				next;
				}

			samples.to.loop <- chips.info$SampleID[this.chip == chips.info$cartridge];
			samples.to.loop <- c(samples.to.loop[!samples.to.loop %in% tmp.ref], tmp.ref);
			}

		# if length of tmp.ref is greater than 1, take the average of the tmp.ref
		if (length(tmp.ref) > 1) {

			# take the average and create a new column in ns.counts as avg.ref
			ns.counts$avg.ref <- apply(X = ns.counts[,tmp.ref], MARGIN = 1, FUN = mean, na.rm = TRUE);

			# and change tmp.ref to 'avg.ref'
			tmp.ref <- 'avg.ref';
			}   

		# avoid divisions by 0 by adding a pseudo-count, if needed
		if (any(na.omit(0 == ns.counts[,tmp.ref]))) { ns.counts[,tmp.ref][0 == ns.counts[,tmp.ref]] <- 1; }

		# divide each test sample probe value by corresponding probes in the ref samples
		output[,samples.to.loop] <- ns.counts[,samples.to.loop] / ns.counts[,tmp.ref];
		}

	# ensure row names are probe names
	rownames(output) <- ns.counts$Name;

	return(output);
	}

