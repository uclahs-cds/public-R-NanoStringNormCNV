collapse.genes <- function(normalized.data){
	# generate normalized.data.aggr gene-level data
	normalized.data.gene <- aggregate(
			x = normalized.data[, !colnames(normalized.data) %in% c('CodeClass', 'Name', 'Accession')],
			by = list(normalized.data$Accession),
			FUN = mean
			);

	# rename the Group.1 column
	colnames(normalized.data.gene)['Group.1' == colnames(normalized.data.gene)] <- 'Name';

	# add gene annotation 
	matching.inds <- unlist(lapply(normalized.data.gene$Name, function(f) which(normalized.data$Accession == f)[1]));
	normalized.data.gene  <- cbind(
		CodeClass = normalized.data[matching.inds, "CodeClass"],
		Name = normalized.data.gene[, "Name"],
		Accession = normalized.data[matching.inds, "Accession"],
		normalized.data.gene[, !colnames(normalized.data.gene) == "Name"]
		);
	rownames(normalized.data.gene) <- normalized.data.gene$Name;

	return(normalized.data.gene);
	}

