
collapse.genes <- function(nano.df){
	# generate nano.df.aggr gene-level data
	nano.df.gene <- aggregate(
			x = nano.df[, !colnames(nano.df) %in% c('CodeClass', 'Name', 'Accession')],
			by = list(nano.df$Accession),
			FUN = mean
			);

	# rename the Group.1 column
	colnames(nano.df.gene)['Group.1' == colnames(nano.df.gene)] <- 'Name';

	return(nano.df.gene);
	}

