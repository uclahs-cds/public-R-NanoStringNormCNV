get.tumour.normal.ratio <- function(normalized.data, reference, chip.info, per.chip = FALSE){
	# create empty data-frame to store data
	output <- normalized.data;

	cols.to.keep   <- colnames(normalized.data);
	cols.to.remove <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	if (any(cols.to.remove %in% cols.to.keep)) {
		cols.to.keep <- cols.to.keep[!cols.to.keep %in% cols.to.remove];
		}
	if (! any(chip.info$SampleID %in% colnames(normalized.data))) {
		chip.info <- chip.info[match(colnames(normalized.data), chip.info$SampleID),];
		}

	output <- output[, cols.to.keep, drop = FALSE];
	output[, cols.to.keep] <- NA;

	# see if user asks to use only chip-specific references
	if (per.chip) {
		chips <- unique(chip.info$Cartridge);
		if (length(chips) < 1) {
			flog.warn("Cannot get tumour/normal ratios per chip: missing cartridge (chip) information!");
			per.chip <- 0;
			}
		}
	if (!per.chip) {
		chips <- 'combined';
		}

	# define sample order here so the reference sample is placed at the end
	samples.to.loop <- c(colnames(output)[!colnames(output) %in% reference], reference);

	# loop over chips
	for (this.chip in chips) {

		# define a temporary reference variable before the start of a new iteration
		tmp.ref <- reference;

		if (this.chip != 'combined') {

			# get the reference(s) for the given chip
			tmp.ref <- chip.info$SampleID[chip.info$Cartridge %in% this.chip & chip.info$SampleID %in% reference];

			# skip when per.chip is requested but there are no reference samples on the chip
			if (length(tmp.ref) < 1) { 
				flog.warn(paste0(
					"No reference samples on cartridge ", this.chip, ".",
					" Try calling CNAs with 'per.chip = FALSE'.",
					" Ratios unavailable for samples:\n\t* ",
					paste(chip.info$SampleID[this.chip == chip.info$Cartridge], collapse = "\n\t* ")
					));
				next;
				}

			samples.to.loop <- chip.info$SampleID[this.chip == chip.info$Cartridge];
			samples.to.loop <- c(samples.to.loop[!samples.to.loop %in% tmp.ref], tmp.ref);
			}

		# if number of reference samples is greater than 1, average the values
		if (length(tmp.ref) > 1) {
			normalized.data$avg.ref <- apply(
				X = normalized.data[,tmp.ref],
				MARGIN = 1,
				FUN = mean,
				na.rm = TRUE
				);
			tmp.ref <- 'avg.ref';
			}   

		# avoid divisions by 0 by adding a pseudo-count where necessary
		if (any(na.omit(0 == normalized.data[,tmp.ref]))) { normalized.data[,tmp.ref][0 == normalized.data[,tmp.ref]] <- 1; }

		# divide each tumour sample probe value by corresponding reference probe values
		output[,samples.to.loop] <- normalized.data[,samples.to.loop] / normalized.data[,tmp.ref];
		}

	# ensure row names are probe names
	rownames(output) <- normalized.data$Name;

	return(output);
	}

