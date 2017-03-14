score.runs <- function(replicate.eval, normalized.data, cna.rounded, phenodata, cna.normals = NULL) {
	# initialize output variable
	scores <- list(
		ari.chip = NA,
		ari.pts.normcor = NA,
		ari.type = NA,
		ari.pts = NA,
		sd.inv = NA,
		sd.hk = NA,
		sd.inv.and.hk = NA,
		replicates.conc = NA,
		prop.disc.genes = NA,
		normals.w.cnas = NA,
		normal.cnas = NA,
		total.cnas = NA
		);

	score.n <- scores;

	# remove unnecessary probes
	normalized.data <- normalized.data[normalized.data$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"),];

	# order annotation
	phenodata <- phenodata[match(colnames(normalized.data[, -(1:3)]), phenodata$SampleID),];

	### Check the adjusted rand index (ARI) of multiple parameters
	# evaluate using cartridge information (all patients)
	flog.info("Calculating ARI using cartridge information (all samples)..");
	scores$ari.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = log10(normalized.data[, -c(1:3)] + 1),
		feature = phenodata$Cartridge,
		is.discrete = FALSE
		);
	score.n$ari.chip <- ncol(normalized.data[, -c(1:3)]);

	# evaluate tumour replicates (normalized NanoString counts)
	if (!is.null(replicate.eval$norm.counts)) {
		replicate.eval$norm.counts <- replicate.eval$norm.counts[,
			replicate.eval$count.pheno$SampleID[which(replicate.eval$count.pheno$Type == 'Tumour')]
			];
		colnames(replicate.eval$norm.counts) <- phenodata[
			match(colnames(replicate.eval$norm.counts), phenodata$SampleID)
			,]$Patient;
		
		flog.info("Calculating ARI using normalized counts (tumour only)..");
		scores$ari.pts.normcor <- NanoStringNormCNV::get.ari(
			data.to.cluster = log10(replicate.eval$norm.counts + 1),
			feature = replicate.eval$count.pheno$Patient[replicate.eval$count.pheno$Type == 'Tumour'],
			is.discrete = FALSE
			);
		score.n$ari.pts.normcor <- ncol(replicate.eval$norm.counts);
		}

	# evaluate using tissue type information (all patients)
	if (any(phenodata$Type == 'Reference')) {
		flog.info("Calculating ARI using tissue type information (all samples)..");
		scores$ari.type <- NanoStringNormCNV::get.ari(
			data.to.cluster = log10(normalized.data[, -c(1:3)] + 1),
			feature = phenodata$Type,
			is.discrete = FALSE
			);
		score.n$ari.type <- ncol(normalized.data[, -c(1:3)]);
		}

	# evaluate tumour replicates (CN calls)
	if (!is.null(replicate.eval$cna.calls)) {
		# must remove probes with NA values from CN calls first
		na.probes <- as.vector(which(is.na(rowSums(replicate.eval$cna.calls))));
		if (length(na.probes) > 0) {
			flog.warn(paste0(
				"Removing the following from 'ari.pts' calculation due to NA values in >=1 samples:\n",
				paste("\t", rownames(replicate.eval$cna.calls)[na.probes], collapse = "\n")
				));
			replicate.eval$cna.calls <- replicate.eval$cna.calls[-na.probes,];
			}

		flog.info("Calculating ARI using copy number calls (tumour only)..");
		scores$ari.pts <- NanoStringNormCNV::get.ari(
			data.to.cluster = replicate.eval$cna.calls,
			feature = phenodata[match(colnames(replicate.eval$cna.calls), phenodata$SampleID),]$Patient
			);
		score.n$ari.pts <- ncol(replicate.eval$cna.calls);
		}

	### Calculate standard deviation for control genes
	scores$sd.inv  <- mean(apply(normalized.data[normalized.data$CodeClass == 'Invariant', -c(1:3)], 1, sd));
	score.n$sd.inv <- ncol(normalized.data[normalized.data$CodeClass == 'Invariant', -c(1:3)]);

	scores$sd.hk  <- mean(apply(normalized.data[normalized.data$CodeClass == 'Housekeeping', -c(1:3)], 1, sd));
	score.n$sd.hk <- ncol(normalized.data[normalized.data$CodeClass == 'Housekeeping', -c(1:3)]);

	scores$sd.inv.and.hk  <- mean(apply(normalized.data[normalized.data$CodeClass %in% c('Housekeeping', 'Invariant'), -c(1:3)], 1, sd));
	score.n$sd.inv.and.hk <- ncol(normalized.data[normalized.data$CodeClass %in% c('Housekeeping', 'Invariant'), -c(1:3)]);

	### Summarize replicate concordance
	if (!is.null(replicate.eval$conc.summary)) {
		scores$replicates.conc <- mean(replicate.eval$conc.summary);
		score.n$replicates.conc <- length(replicate.eval$conc.summary);
		}
	
	if (!is.null(replicate.eval$concordance)) {
		scores$prop.disc.genes <- sum(
			as.vector(apply(
				X = replicate.eval$concordance, 
				MARGIN = 1, 
				FUN = function(f) ifelse(any(f == 0), 1, 0)
				))
			)/nrow(replicate.eval$concordance);
		score.n$prop.disc.genes <- ncol(replicate.eval$concordance);
		}

	### Summarize CNA counts
	if (!is.null(cna.normals)) {
		# number of normal samples where at least one CNA was called
		scores$normals.w.cnas  <- length(which(apply(cna.normals, 2, function(f) any(f != 2))));
		score.n$normals.w.cnas <- ncol(cna.normals);
		# total number of CNAs called across all normal samples
		scores$normal.cnas <-  sum(apply(cna.normals, 1, function(f)  sum(f != 2, na.rm = TRUE)), na.rm = TRUE);
		score.n$normal.cnas <- ncol(cna.normals);
		}

	# total number of CNAs called across all tumour samples
	scores$total.cnas  <- sum(apply(cna.rounded, 2, function(f) sum(f != 2, na.rm = TRUE)), na.rm = TRUE);
	score.n$total.cnas <- ncol(cna.rounded);

	return(list(scores = scores, sample.size = score.n));
	}
