normalize.per.chip <- function(phenodata, raw.data, cc, bc, sc, oth, do.nsn, do.rcc.inv, covs, transform.data = TRUE, plot.types = 'all'){
	# modify header to NanoStringNorm standard
	colnames(phenodata)[colnames(phenodata) == 'Cartridge'] <- 'cartridge';
	colnames(phenodata)[colnames(phenodata) == 'Type'] 		<- 'type';

	nano.parts <- list();
	cartridges <- unique(phenodata$cartridge);

	rownames(raw.data) <- raw.data$Name;

	for (chip in 1:length(cartridges)) {
		flog.info(paste("Normalizing cartridge", cartridges[chip]), "probes..");
		cur.samples <- which(phenodata$cartridge == cartridges[chip]);

		if (length(unique(covs[cur.samples, 'type'])) > 1) {
			use.covs <- covs[cur.samples, 'type', drop = FALSE];
		} else {
			use.covs  <- NA;
			}

		# normalization using code count, background noise, sample content
		if (do.nsn) {
			nano.parts[[chip]] <- NanoStringNorm::NanoStringNorm(
				x = raw.data[,(cur.samples + 3), drop = FALSE],
				CodeCount = cc,
				Background = bc,
				SampleContent = sc,
				round.values = FALSE,
				take.log = FALSE,
				traits = use.covs,
				anno = raw.data[, 1:3]
				);
			
			pdf(paste0('NanoStringNorm_plots_chip', chip, '.pdf'));
			NanoStringNorm::Plot.NanoStringNorm(
				x = nano.parts[[chip]],
				label.best.guess = TRUE,
				plot.type = plot.types
				);
			dev.off();
			
			nano.parts[[chip]] <- nano.parts[[chip]]$normalized.data;
			colnames(nano.parts[[chip]])[1] <- 'CodeClass';
		} else {
			nano.parts[[chip]] <- raw.data[,c(1:3, (cur.samples + 3)), drop = FALSE];
			}

		# invariant probe normalization
		if (do.rcc.inv) {
			if (all(c('SampleID', 'type') %in% colnames(phenodata))) {
				phenodata.inv <- phenodata[,c('SampleID', 'type')];
			} else {
				phenodata.inv <- NULL;
				}
			colnames(phenodata.inv)[colnames(phenodata.inv) == 'type'] <- 'Type';
			nano.parts[[chip]] <- NanoStringNormCNV::invariant.probe.norm(nano.parts[[chip]], phenodata.inv);
			}

		# perform 'other' normalization last
		if (do.nsn & oth != 'none') {
			nano.parts[[chip]] <- NanoStringNorm::NanoStringNorm(
				x = nano.parts[[chip]][, -(1:3)],
				OtherNorm = oth,
				round.values = FALSE,
				take.log = FALSE,
				traits = use.covs,
				anno = raw.data[, 1:3]
				);
			nano.parts[[chip]] <- nano.parts[[chip]]$normalized.data;
			colnames(nano.parts[[chip]])[1] <- 'CodeClass';
			}
		}

	# combine data
	normalized.data <- cbind(raw.data[, 1:3], do.call(cbind, lapply(nano.parts, function(f) f[, -c(1:3), drop = FALSE])));
	normalized.data <- normalized.data[, c(1:3, order(colnames(normalized.data)[-(1:3)]) + 3)];
	rownames(normalized.data) <- raw.data[, colnames(raw.data) == 'Name'];

	# transform data so there are no negative counts
	if (transform.data & any(normalized.data[, -(1:3)] < 0)) {
		normalized.data[, -(1:3)] <- (normalized.data[, -(1:3)] - min(normalized.data[, -(1:3)])) + 0.1;
		}

	return(normalized.data);
	}
