collapse.genes <- function(normalized.data){
	# split control and non-control probes
	normalized.data.non.ctrl <- normalized.data[ normalized.data$CodeClass %in% c('Endogenous', 'Housekeeping'),];
	normalized.data.ctrl     <- normalized.data[!normalized.data$CodeClass %in% c('Endogenous', 'Housekeeping'),];

	# generate normalized.data.aggr gene-level data for non-control probes
	normalized.data.gene <- aggregate(
		x = normalized.data.non.ctrl[, !colnames(normalized.data.non.ctrl) %in% c('CodeClass', 'Name', 'Accession')],
		by = list(normalized.data.non.ctrl$Accession),
		FUN = mean
		);

	# rename the Group.1 column
	colnames(normalized.data.gene)['Group.1' == colnames(normalized.data.gene)] <- 'Name';

	# add gene annotation 
	matching.inds <- unlist(lapply(normalized.data.gene$Name, function(f) which(normalized.data.non.ctrl$Accession == f)[1]));
	normalized.data.gene  <- cbind(
		CodeClass = normalized.data.non.ctrl[matching.inds, "CodeClass"],
		Name = normalized.data.gene[, "Name"],
		Accession = normalized.data.non.ctrl[matching.inds, "Accession"],
		normalized.data.gene[, !colnames(normalized.data.gene) == "Name"],
		stringsAsFactors = FALSE
		);
	rownames(normalized.data.gene) <- normalized.data.gene$Name;

	# re-combine with (non-aggregated) control probes
	normalized.data.gene <- rbind(normalized.data.gene, normalized.data.ctrl);

	return(normalized.data.gene);
	}

