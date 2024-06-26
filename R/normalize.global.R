normalize.global <- function(raw.data, cc, bc, sc, oth, do.rcc.inv, covs, transform.data = TRUE, plot.types = 'all', phenodata = NULL){
	# modify header to NanoStringNorm standard
	colnames(phenodata)[colnames(phenodata) == 'Cartridge'] <- 'cartridge';
	colnames(phenodata)[colnames(phenodata) == 'Type'] 		<- 'type';

	# normalization using code count, background noise, sample content
	if (cc != 'none' | bc != 'none' | sc != 'none') {
		nano.norm <- NanoStringNorm::NanoStringNorm(
			x = raw.data[, -c(1:3)],
			CodeCount = cc,
			Background = bc,
			SampleContent = sc,
			round.values = FALSE,
			take.log = FALSE,
			traits = covs,
			anno = raw.data[, 1:3]
			);
		normalized.data <- nano.norm$normalized.data;
		colnames(normalized.data)[1] <- 'CodeClass';

		pdf('NanoStringNorm_plots_all.pdf');
		NanoStringNorm::Plot.NanoStringNorm(
			x = nano.norm,
			label.best.guess = TRUE,
			plot.type = plot.types
			);
		dev.off();
	} else {
		normalized.data <- raw.data;
		}

	# invariant probe normalization
	if (do.rcc.inv) {
		colnames(phenodata)[colnames(phenodata) == 'type'] <- 'Type';
		normalized.data <- NanoStringNormCNV::invariant.probe.norm(normalized.data, phenodata);
		}

	# perform 'other' normalization last
	if (oth != 'none') {
		nano.norm <- NanoStringNorm::NanoStringNorm(
			x = normalized.data[, -(1:3)],
			OtherNorm = oth,
			round.values = FALSE,
			take.log = FALSE,
			traits = covs,
			anno = normalized.data[, 1:3]
			);
		normalized.data <- nano.norm$normalized.data;
		colnames(normalized.data)[1] <- 'CodeClass';
		}

	# transform data so there are no negative counts
	if (transform.data & any(normalized.data[, -(1:3)] < 0)) {
		normalized.data[, -(1:3)] <- (normalized.data[, -(1:3)] - min(normalized.data[, -(1:3)])) + 0.1;
		}

	# ensure row names are probe names
	rownames(normalized.data) <- raw.data[, colnames(raw.data) == 'Name'];

	return(normalized.data);
	}
