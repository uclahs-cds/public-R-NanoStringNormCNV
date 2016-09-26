
collapse.duplicate.samples <- function(nano.df, annot){
	# for each sample that are duplicated, merge the counts so you get a single count
	nano.df.aggr <- nano.df[,c('CodeClass', 'Name', 'Accession')];

	# loop over each sample
	for (sample.name in unique(annot$Name)) {

		# get the colnames equivalent of sample.name
		sample.id <- annot$SampleID[sample.name == annot$Name];

		# aggregate this by median value
		nano.df.aggr[,sample.name] <- apply(
			X = nano.df[,sample.id],
			MARGIN = 1,
			FUN = median
			);
		}

	return(nano.df.aggr);
	}

