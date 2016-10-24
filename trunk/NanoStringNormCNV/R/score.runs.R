score.runs <- function(replicates, normalized, cnas, sample.annot, normals = NULL) {
	scores <- list();

	# remove unnecessary probes
	normalized <- normalized[normalized$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"),];

	# order annotation
	sample.annot <- sample.annot[match(colnames(normalized[, -(1:3)]), sample.annot$SampleID),];

	### Check the adjusted rand index (ARI) of multiple parameters
	# evaluate batches: make the correlation matrix using all samples and, again, for tumours only
	scores$ari.chip <- NanoStringNormCNV::get.ari(
		data = log10(normalized[, -c(1:3)] + 1),
		feature = sample.annot$cartridge
		);

	replicates$norm.counts <- replicates$norm.counts[, replicates$count.pheno$SampleID[which(replicates$count.pheno$type == 'Tumour')]];
	colnames(replicates$norm.counts) <- phenodata[match(colnames(replicates$norm.counts), phenodata$SampleID),]$Patient;
	
	scores$ari.pts.normcor <- NanoStringNormCNV::get.ari(
		data = log10(replicates$norm.counts + 1),
		feature = replicates$count.pheno$Patient[replicates$count.pheno$type == 'Tumour']
		);

	# type, use all patients
	if (any(sample.annot$type == 'Reference')) {
		scores$ari.type <- NanoStringNormCNV::get.ari(
			data = log10(normalized[, -c(1:3)] + 1),
			feature = sample.annot$type
			);
		}

	# must remove CNA call probes that contain any NAs
	na.probes <- as.vector(which(is.na(rowSums(replicates$cna.calls))));
	if (length(na.probes) > 0) {
		flog.warn(paste0(
			"Removing the following from 'ari.pts' calculation due to NA values in one or more samples:\n",
			paste("\t", rownames(replicates$cna.calls)[na.probes], collapse = "\n")
			));
		replicates$cna.calls <- replicates$cna.calls[-na.probes,];
		}

	# replicates using CNA data
	scores$ari.pts <- NanoStringNormCNV::get.ari(
		data = replicates$cna.calls,
		feature = phenodata[match(colnames(replicates$cna.calls), phenodata$SampleID),]$Patient
		);

	### SD for control genes
	scores$sd.inv 		 <- mean(apply(normalized[normalized$CodeClass == 'Invariant', -c(1:3)], 1, sd));
	scores$sd.hk  		 <- mean(apply(normalized[normalized$CodeClass == 'Housekeeping', -c(1:3)], 1, sd));
	scores$sd.inv.and.hk <- mean(apply(normalized[normalized$CodeClass == 'Housekeeping' | normalized$CodeClass == 'Invariant', -c(1:3)], 1, sd));

	### Replicate concordance
	scores$replicates.conc <- mean(replicates$conc.summary);
	scores$prop.disc.genes <- sum(as.vector(apply(replicates$concordance, 1, function(f) ifelse(any(f == 0), 1, 0))))/nrow(replicates$concordance);

	# check how many CNAs were called in the normal samples
	if (!is.null(normals)) {
		scores$normals.w.cnas <- length(which(apply(normals, 2, function(f) any(f != 2))));
		scores$normal.cnas 	  <- sum(apply(normals, 1, function(f)  sum(f != 2, na.rm = TRUE)), na.rm = TRUE);
	} else {
		scores$normals.w.cnas <- NA;
		scores$normal.cnas <- NA;
		}

	# check how many CNAs were called in the tumour samples
	scores$total.cnas <- sum(apply(cnas, 2, function(f) sum(f != 2, na.rm = TRUE)), na.rm = TRUE);

	return(scores);
	}
