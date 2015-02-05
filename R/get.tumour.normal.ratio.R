
get.tumour.normal.ratio <- function(ref, nchips, chips.annot, ns.counts, output){
	# define samples.to.loop here so the reference sample is placed at the end
	samples.to.loop <- c(colnames(output)[!colnames(output) %in% ref], ref);

	# loop over nchips
	for (perm.i in 1:nchips) {

		# define a tmp.ref every time before the start of a new iteration
		tmp.ref <- ref;

		# if per.chip is TRUE
		if (per.chip) {

			# get the tmp.ref for chip specifically
			this.chip <- paste0('Chip ', perm.i);
			tmp.ref <- chips.annot$SampleID[chips.annot$Chip %in% this.chip & chips.annot$SampleID %in% ref];

			# when per chip is requested, but there are no ref samples on the chip, use pooled refs then
			if (length(tmp.ref) < 1) { next; }

			samples.to.loop <- chips.annot$SampleID[this.chip == chips.annot$Chip];
			samples.to.loop <- c(samples.to.loop[!samples.to.loop %in% tmp.ref], tmp.ref);
			}

		# if length of tmp.ref is greater than 1, take the average of the tmp.ref
		if (length(tmp.ref) > 1) {

			# take the average and create a new column in ns.counts as avg.ref
			ns.counts$avg.ref <- apply(X = ns.counts[,tmp.ref], MARGIN = 1, FUN = mean);

			# and change tmp.ref to 'avg.ref'
			tmp.ref <- 'avg.ref';
			}

		# loop over each sample first and then do the processing for reference samples!
		for (this.sample in samples.to.loop) {

			# divide each test samples probe value by corresponding probes in the ref samples
			if (any(0 == ns.counts[,tmp.ref])) { ns.counts[,tmp.ref][0 == ns.counts[,tmp.ref]] <- 1; }
			tmp.ratio <- ns.counts[,this.sample] / ns.counts[,tmp.ref];
			output[,this.sample] <- tmp.ratio;
			}
		}

	return(output);
	}

