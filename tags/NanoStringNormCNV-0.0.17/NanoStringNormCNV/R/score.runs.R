
score.runs <- function(replicates, normalized, cnas, sample.annot, genes, collapsed, normals) {
	scores <- list();

	# determine if the normalized data is counts or continuous (actually can use > or < 0)
	discrete.data <- FALSE;
	if (min(normalized[, -c(1:3)]) >= 0) {
		discrete.data <- TRUE;
		}

	### Check the adjusted rand index (ARI) of multiple parameters
	
	# evaluate batches
	# make the correlation matrix using all samples
	# and again for tumours only
	if (discrete.data) {
		scores$ari.chip <- NanoStringNormCNV::get.ari(
			log10(normalized[, -c(1:3)] + 1),
			sample.annot$cartridge
			);

		replicates$norm.counts <- replicates$norm.counts[ , replicates$pheno$SampleID[which(replicates$pheno$type == 'Tumour')]];
		colnames(replicates$norm.counts) <- substr(colnames(replicates$norm.counts), 1, 8);
		
		scores$ari.pts.normcor <- NanoStringNormCNV::get.ari(
			log10(replicates$norm.counts+1),
			replicates$pheno$Patient[replicates$pheno$type == 'Tumour']
			);
	} else {
		scores$ari.chip <- NanoStringNormCNV::get.ari(
			normalized[, -c(1:3)] + 1,
			sample.annot$cartridge
			);

		replicates$norm.counts <- replicates$norm.counts[ , replicates$pheno$SampleID[which(replicates$pheno$type == 'Tumour')]];
		colnames(replicates$norm.counts) <- substr(colnames(replicates$norm.counts), 1, 8);
		
		scores$ari.pts.normcor <- NanoStringNormCNV::get.ari(
			replicates$norm.counts,
			replicates$pheno$Patient[replicates$pheno$type == 'Tumour']
			);
		}

	# type, use all patients
	if (any(sample.annot$type == 'Reference')) {
		if (discrete.data) {
			scores$ari.type <- NanoStringNormCNV::get.ari(log10(normalized[, -c(1:3)] + 1), sample.annot$type);
		} else {
			scores$ari.type <- NanoStringNormCNV::get.ari(normalized[, -c(1:3)] + 1, sample.annot$type);
			}
		}

	# replicates using CNA data
	scores$ari.pts <- NanoStringNormCNV::get.ari(
		t(replicates$cnas),
		substr(colnames(replicates$cnas), 1, 8),
		discrete.data = TRUE
		);

	### SD for control genes
	if (discrete.data) {
		scores$sd.inv 		 <- mean(apply(normalized[normalized$CodeClass == 'Invariant', -c(1:3)], 1, sd));
		scores$sd.hk  		 <- mean(apply(normalized[normalized$CodeClass == 'Housekeeping', -c(1:3)], 1, sd));
		scores$sd.inv.and.hk <- mean(apply(normalized[normalized$CodeClass == 'Housekeeping' | normalized$CodeClass == 'Invariant', -c(1:3)], 1, sd));
		}

	### Replicate concordance
	scores$replicates.conc <- mean(replicates$conc.summary);
	scores$prop.disc.genes <- sum(as.vector(apply(replicates$concordance, 1, function(f) ifelse(any(f == 0), 1, 0))))/nrow(replicates$concordance);

	# check how many CNAs were called in the normal samples
	print(head(normals));
	scores$normals.w.cnas <- length(which(apply(normals, 2, function(f) any(f !=2 ))));
	scores$normal.cnas 	  <- sum(apply(normals, 1, function(f)  sum(f != 2)));

	# check how many CNAs were called in the tumour samples
	scores$total.cnas <- sum(apply(cnas, 2, function(f) sum(f !=2 )));

	return(scores);
	}
