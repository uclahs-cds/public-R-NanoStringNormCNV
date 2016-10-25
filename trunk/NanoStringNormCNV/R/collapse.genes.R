collapse.genes <- function(nano.df){
	# generate nano.df.aggr gene-level data
	nano.df.gene <- aggregate(
			x = nano.df[, !colnames(nano.df) %in% c('CodeClass', 'Name', 'Accession')],
			by = list(nano.df$Accession),
			FUN = mean
			);

	# rename the Group.1 column
	colnames(nano.df.gene)['Group.1' == colnames(nano.df.gene)] <- 'Name';

	# add gene annotation 
	matching.inds <- unlist(lapply(nano.df.gene$Name, function(f) which(nano.df$Accession == f)[1]));
	nano.df.gene  <- cbind(
		CodeClass = nano.df[matching.inds, "CodeClass"],
		Name = nano.df.gene[, "Name"],
		Accession = nano.df[matching.inds, "Accession"],
		nano.df.gene[, !colnames(nano.df.gene) == "Name"]
		);
	rownames(nano.df.gene) <- nano.df.gene$Name;

	return(nano.df.gene);
	}

