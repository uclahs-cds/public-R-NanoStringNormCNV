
calculate.replicate.concordance <- function (input, pheno, gene.names) {

	# initiate an out.table file
	out.table <- as.data.frame(
		matrix(
			nrow = nrow(input),
			ncol = length(unique(pheno$Name)),
			dimnames = list(
				gene.names,
				unique(pheno$Name)
				)
			)
		);

	# loop over each unique sample from Name and create a difference matrix
	for (this.sample in unique(pheno$Name)) {
		# get the samples pertaining to the name
		this.ID <- pheno[pheno$Name == this.sample, 'SampleID'];

		# calculate the per-gene variance of replicates
		per.gene.conc <- apply(
			X = input[,which(colnames(input) == this.ID), drop = FALSE],
			MARGIN = 1,
			FUN = function(f)
				ifelse(length(unique(as.numeric(f))) == 1, 1, 0)
			);

		# add that to the out.table
		out.table[,this.sample] <- per.gene.conc;
		}

	# return out.table
	return (out.table);
	}

