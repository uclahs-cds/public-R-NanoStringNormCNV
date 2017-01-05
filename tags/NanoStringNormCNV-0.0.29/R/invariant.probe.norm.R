invariant.probe.norm <- function(normalized.data, phenodata = NULL){
	# NB: it is not necessary to normalize the data for POS, NEG or Restriction Site controls
	# INV-normalization procedure should be carried out on the Inv and End code classes.

	# remove CodeClass, Name, Accession from colnames
	grep.cols <- colnames(normalized.data)[!colnames(normalized.data) %in% c('CodeClass', 'Name', 'Accession')];

	# check if sample ID and tissue type info exists
	if (!is.null(phenodata)) {
		if (all(c('SampleID', 'Type') %in% colnames(phenodata))) { phenodata <- phenodata[,c('SampleID', 'Type')]; }
	} else {
		phenodata <- NULL;
		}

	# step 1
	# get average values for all the invariant probes 
	inv.avg <- apply(
		X = normalized.data['Invariant' == normalized.data$CodeClass, grep.cols, drop = FALSE],
		MARGIN = 2,
		FUN = mean
		);

	# step 2A
	# calculate mean of average INV count values
	mean.inv.avg <- mean(inv.avg);

	# step 2B
	# plot the counts for the invariant probes. It is BAD to have counts < 100 for these probes as it signifies
	# low DNA input, leading to unreliable CNA calls. Especially if the low counts are in the reference samples!!
	NanoStringNormCNV::make.invariant.probe.plot(
		inv.probe.counts = normalized.data['Invariant' == normalized.data$CodeClass, grep.cols, drop = FALSE],
		tissue.type = phenodata
		);

	# step 3
	# calculate normalization factor
	norm.factor <- mean.inv.avg / inv.avg;

	# create a new object for nano.norm
	rows.to.norm <- 1:nrow(normalized.data);
	nano.norm    <- normalized.data[rows.to.norm, c('CodeClass', 'Name', 'Accession')];
	nano.norm[,grep.cols] <- NA;

	# sanity check to make sure all names are the same
	if (!all(names(norm.factor) == grep.cols)) { stop ('Mismatch in Names!'); }

	# loop over each sample
	for (ind in 1:length(grep.cols)) {

		# get the ind-th row
		this.col    <- normalized.data[rows.to.norm, grep.cols[ind], drop = FALSE];
		this.factor <- norm.factor[ind];

		# get INV-normalized counts (RAW counts for each probe x normalization factor)
		nano.norm[,grep.cols[ind]] <- this.col * this.factor;
		}

	return(nano.norm);
	}
