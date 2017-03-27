invariant.probe.norm <- function(nanostring.data, phenodata = NULL){
	# remove CodeClass, Name, Accession from colnames
	grep.cols <- colnames(nanostring.data)[!colnames(nanostring.data) %in% c('CodeClass', 'Name', 'Accession')];

	# check if sample ID and tissue type info exists
	if (!is.null(phenodata) & !all(c('SampleID', 'Type') %in% colnames(phenodata))) {
		phenodata <- NULL;
		}

	# get average values for all the invariant probes, per sample
	inv.avg <- apply(
		X = nanostring.data['Invariant' == nanostring.data$CodeClass, grep.cols, drop = FALSE],
		MARGIN = 2,
		FUN = mean
		);

	# calculate mean of average INV count values
	mean.inv.avg <- mean(inv.avg);

	# plot counts of invariant probes: counts of < 100 signify low DNA input, leading to unreliable CNA calls
	# this is especially problematic if the low counts are in the reference samples!
	NanoStringNormCNV::make.invariant.probe.plot(
		inv.probe.counts = nanostring.data['Invariant' == nanostring.data$CodeClass, grep.cols, drop = FALSE],
		tissue.type = phenodata
		);

	# calculate normalization factor
	norm.factor <- mean.inv.avg / inv.avg;

	# create a new object where 'Positive', 'Negative' and 'RestrictionSite' probes will not be normalized
	rows.to.norm <- !(nanostring.data$CodeClass %in% c('Positive', 'Negative', 'RestrictionSite'));
	nano.norm    <- nanostring.data[, c('CodeClass', 'Name', 'Accession', grep.cols)];
	nano.norm[rows.to.norm, grep.cols] <- NA;

	# sanity check to make sure all names are the same
	if (!all(names(norm.factor) == grep.cols)) { stop('Mismatch in Names!'); }

	# get INV-normalized counts (counts for each probe x normalization factor) per sample
	for (ind in 1:length(grep.cols)) {
		this.col    <- nanostring.data[rows.to.norm, grep.cols[ind], drop = FALSE];
		this.factor <- norm.factor[ind];

		nano.norm[rows.to.norm, grep.cols[ind]] <- this.col * this.factor;
		}

	return(nano.norm);
	}
