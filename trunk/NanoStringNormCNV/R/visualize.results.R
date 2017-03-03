visualize.results <- function(raw.data, normalized.data, phenodata = NULL, cna.rounded = NULL, cna.raw = NULL, max.cn = 5, replicate.eval = NULL, exclude.covs = FALSE) {
	# set up covariates
	sample.covs <- gene.covs.pre.norm <- gene.covs.post.norm <- NULL;

	# process covariate information
	if (!exclude.covs) {
		if (! is.null(phenodata)) {
			sample.covs <- phenodata[, colnames(phenodata) %in% c('SampleID', 'Type', 'Cartridge'), drop = FALSE];
			if (! any(colnames(sample.covs) %in% 'SampleID')) {
				stop("Must include sample IDs in phenodata!");
				}
			if (ncol(sample.covs) < 2) {
				sample.covs <- NULL;
				flog.warn("Excluding sample covariates due to missing information!");
				}
			}
		gene.covs.pre.norm  <- raw.data[, c('Name', 'CodeClass')];
		gene.covs.post.norm <- normalized.data[, c('Name', 'CodeClass')];
		}

	# re-format plotting data
	cols.to.remove <- c('CodeClass', 'Name', 'Accession');
	
	raw.data.reformatted  <- raw.data[,  !(colnames(raw.data)  %in% cols.to.remove)];
	normalized.data.reformatted <- normalized.data[, !(colnames(normalized.data) %in% cols.to.remove)];
	
	row.names(raw.data.reformatted)  <- raw.data$Name;
	row.names(normalized.data.reformatted) <- normalized.data$Name;

	if (!is.null(cna.raw)) { cna.raw[cna.raw > max.cn] <- max.cn; }

	replicate.eval$count.pheno$Patient <- factor(replicate.eval$count.pheno$Patient);
	replicate.eval$count.pheno$Type    <- factor(replicate.eval$count.pheno$Type);
	if (any(names(replicate.eval$count.pheno) == 'outlier')) {
		replicate.eval$count.pheno$outlier <- factor(replicate.eval$count.pheno$outlier, levels = c(0,1));
		}

	replicate.eval$raw.counts <- raw.data[, replicate.eval$count.pheno$SampleID];
	rownames(replicate.eval$raw.counts) <- raw.data$Name;

	##############################
	### Heatmaps (all samples) ###
	##############################

	# raw NanoString counts
	if (nlevels(as.factor(unlist(raw.data.reformatted))) > 1) {
		flog.info("Plotting raw counts heatmap..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = raw.data.reformatted,
			fname.stem = 'raw',
			covs.rows = sample.covs,
			covs.cols = gene.covs.pre.norm
			);
	} else {
		flog.warn("Unable to plot raw counts data!");
		}

	# normalized NanoString counts
	if (nlevels(as.factor(unlist(normalized.data.reformatted))) > 1) {
		flog.info("Plotting normalized counts heatmap..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = normalized.data.reformatted,
			fname.stem = 'normalized',
			covs.rows = sample.covs,
			covs.cols = gene.covs.post.norm
			);
	} else {
		flog.warn("Unable to plot normalized count data!");
		}

	# raw count correlations
	if (nlevels(as.factor(unlist(raw.data.reformatted))) > 1) {
		flog.info("Plotting raw count correlations heatmap..");
		NanoStringNormCNV::make.sample.correlations.heatmap(
			nano.counts = log10(raw.data.reformatted + 1),
			fname.stem = 'raw-count',
			covs = sample.covs
			);
	} else {
		flog.warn("Unable to plot raw count correlations!");
		}

	# normalized count correlations
	if (nlevels(as.factor(unlist(normalized.data.reformatted))) > 1) {
		flog.info("Plotting normalized count correlations heatmap..");
		NanoStringNormCNV::make.sample.correlations.heatmap(
			nano.counts = log10(normalized.data.reformatted + 1),
			fname.stem = 'norm-count',
			covs = sample.covs
			);
	} else {
		flog.warn("Unable to plot normalized count correlations!");
		}

	# rounded CNA calls
	if (! is.null(cna.rounded)){
		if (nlevels(as.factor(unlist(cna.rounded))) > 1) {
			flog.info("Plotting 'rounded CNA calls' heatmap..");
			NanoStringNormCNV::make.cna.heatmap(
				nano.cnas = cna.rounded,
				fname.stem = 'rounded-cna-calls',
				rounded = TRUE,
				covs.cols = gene.covs.post.norm,
				covs.rows = sample.covs,
				min.cn = 0,
				width = 7
				);
		} else {
			flog.warn("Unable to plot 'rounded CNA' call data!");
			}
		}

	# raw CNA calls
	if (! is.null(cna.raw)) {
		if (nlevels(as.factor(unlist(cna.raw))) > 1) {
			flog.info("Plotting 'raw CNA calls' heatmap..");
			NanoStringNormCNV::make.cna.heatmap(
				nano.cnas = cna.raw,
				fname.stem = 'raw-cna-calls',
				rounded = TRUE,
				covs.cols = gene.covs.post.norm,
				covs.rows = sample.covs,
				min.cn = 0,
				width = 7
				);
		} else {
			flog.warn("Unable to plot 'raw CNA' call data!");
			}
		}

	###################################
	### Density plots (all samples) ###
	###################################
	
	# rounded CNA calls
	if (! is.null(cna.rounded)) {
		flog.info("Plotting 'rounded CNA calls' density plots..");
		NanoStringNormCNV::make.cna.densities.plots(
			nano.cnas = cna.rounded,
			fname.stem = 'rounded-cna-calls'
			);
		}

	# raw CNA calls
	if (! is.null(cna.raw)) {
		flog.info("Plotting 'raw CNA calls' density plots..");
		NanoStringNormCNV::make.cna.densities.plots(
			nano.cnas = cna.raw,
			fname.stem = 'raw-cna-calls'
			);
		}

	##################################
	### Heatmaps (replicates only) ###
	##################################

	if (! is.null(replicate.eval)) {
		# raw NanoString counts
		if (nlevels(as.factor(unlist(replicate.eval$raw.counts))) > 1) {
			flog.info("Plotting raw counts heatmap for replicates..");
			NanoStringNormCNV::make.counts.heatmap(
				nano.counts = replicate.eval$raw.counts,
				fname.stem = 'replicate_raw',
				covs.rows = sample.covs,
				covs.cols = gene.covs.pre.norm,
				);
		} else {
			flog.warn("Unable to plot raw count data for replicates!");
			}
		
		# normalized NanoString counts
		if (nlevels(as.factor(unlist(replicate.eval$norm.counts))) > 1) {
			flog.info("Plotting normalized counts heatmap for replicates..");
			NanoStringNormCNV::make.counts.heatmap(
				nano.counts = replicate.eval$norm.counts,
				fname.stem = 'replicate_norm',
				covs.cols = gene.covs.post.norm,
				covs.rows = sample.covs
				);
		} else {
			flog.warn("Unable to plot normalized count data for replicates!");
			}

		# CNA calls
		if (nlevels(as.factor(unlist(replicate.eval$cna.calls))) > 1) {
			flog.info("Plotting CNA calls heatmap for replicates..");
			NanoStringNormCNV::make.cna.heatmap(
				nano.cnas = replicate.eval$cna.calls,
				fname.stem = 'replicate_cna-calls',
				rounded = TRUE,
				covs.cols = gene.covs.post.norm,
				covs.rows = sample.covs,
				min.cn = 0
				);
		} else {
			flog.warn("Unable to plot CNA calls for replicates!");
			}


		# replicate CNA call concordance
		if (nlevels(as.factor(unlist(replicate.eval$concordance))) > 1) {
			flog.info("Plotting replicate CNA call concordance heatmap..");
			NanoStringNormCNV::make.counts.heatmap(
				nano.counts = replicate.eval$concordance,
				fname.stem = 'replicate_cna-concordance',
				covs.cols = gene.covs.post.norm,
				print.ylab = TRUE
				);
		} else {
			flog.warn("Unable to plot CNA call concordance for replicates!");
			}

		# normalized count correlations
		if (nlevels(as.factor(unlist(replicate.eval$norm.counts))) > 1) {
			flog.info("Plotting normalized count correlations heatmap for replicates..");
			NanoStringNormCNV::make.sample.correlations.heatmap(
				nano.counts = log10(replicate.eval$norm.counts + 1),
				fname.stem = 'replicate_norm-count',
				covs = sample.covs
				);
		} else {
			flog.warn("Unable to plot normalized count correlations for replicates!");
			}

		# normalized count correlations (tumour only)
		norm.counts.tumour.only <- replicate.eval$norm.counts[, which(replicate.eval$count.pheno$Type == 'Tumour')];
		if (nlevels(as.factor(unlist(norm.counts.tumour.only))) > 1) {
			flog.info("Plotting normalized count correlations heatmap for replicates (tumour only)..");
			NanoStringNormCNV::make.sample.correlations.heatmap(
				nano.counts = log10(norm.counts.tumour.only + 1),
				fname.stem = 'replicate_tumour-only_norm-count',
				covs = sample.covs
				);
		} else {
			flog.warn("Unable to plot normalized count correlations for tumour only replicates!");
			}
		}
	}
