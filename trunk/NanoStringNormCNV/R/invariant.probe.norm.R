
invariant.probe.norm <- function(nano.df, pheno = NULL){

	# NB: it is not necessary to normalize the data for POS, NEG or Restriction Site controls
	# INV-normalization procedure should be carried out on the Inv and End code classes.

	# remove CodeClass, Name, Accession from colnames
	grep.cols <- colnames(nano.df)[!colnames(nano.df) %in% c('CodeClass', 'Name', 'Accession')];

	# check if sample ID and type info exists
	if (!is.null(pheno)) {
		if (all(c('SampleID', 'type') %in% colnames(pheno))) { pheno <- pheno[,c('SampleID', 'type')]; }
	} else {
		pheno <- NULL;
		}

	#--- invariant probe normalization ---------------------------------------------------------------#
	# step 1
	# get average values for all the invariant probes 
	inv.avg <- apply(
		X = nano.df['Invariant' == nano.df$CodeClass, grep.cols, drop = FALSE],
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
		nano.df['Invariant' == nano.df$CodeClass, grep.cols, drop = FALSE],
		fname.stem = 'invariant_probe',
		sample.type = pheno
		);

	# step 3
	# calculate normalization factor by mean.inv.avg / inv.avg
	norm.factor <- mean.inv.avg / inv.avg;

	# create a new object for nano.norm
	rows.to.norm <- 1:nrow(nano.df);
	nano.norm    <- nano.df[rows.to.norm, c('CodeClass', 'Name', 'Accession')];
	nano.norm[,grep.cols] <- NA;

	# sanity check to make sure all names are the same
	if (!all(names(norm.factor) == grep.cols)) { stop ('Mismatch in Names!'); }

	# loop over each sample
	for (ind in 1:length(grep.cols)) {

		# get the ind-th row
		this.col    <- nano.df[rows.to.norm, grep.cols[ind], drop = FALSE];
		this.factor <- norm.factor[ind];

		# get INV-normalized counts (RAW counts for each probe x normalization factor)
		nano.norm[,grep.cols[ind]] <- this.col * this.factor;
		}

	return(nano.norm);
	}
