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
	flog.info("Plotting raw counts heatmap..");
	NanoStringNormCNV::make.counts.heatmap(
		nano.counts = raw.data.reformatted,
		fname.stem = 'raw',
		covs.rows = sample.covs,
		covs.cols = gene.covs.pre.norm
		);

	# normalized NanoString counts
	flog.info("Plotting normalized counts heatmap..");
	NanoStringNormCNV::make.counts.heatmap(
		nano.counts = normalized.data.reformatted,
		fname.stem = 'normalized',
		covs.rows = sample.covs,
		covs.cols = gene.covs.post.norm
		);

	# raw count correlations
	flog.info("Plotting raw count correlations heatmap..");
	NanoStringNormCNV::make.sample.correlations.heatmap(
		nano.counts = log10(raw.data.reformatted + 1),
		fname.stem = 'raw-count',
		covs = sample.covs
		);

	# normalized count correlations
	flog.info("Plotting normalized count correlations heatmap..");
	NanoStringNormCNV::make.sample.correlations.heatmap(
		nano.counts = log10(normalized.data.reformatted + 1),
		fname.stem = 'norm-count',
		covs = sample.covs
		);

	# rounded CNA calls
	if (! is.null(cna.rounded)){
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
		}

	# raw CNA calls
	if (! is.null(cna.raw)) {
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
		}

	###################################
	### Density plots (all.samples) ###
	###################################

	# raw CNA calls
	if (! is.null(cna.raw)) {
		flog.info("Plotting 'raw CNA calls' density plots..");
		NanoStringNormCNV::make.cna.densities.plots(
			nano.cnas = cna.raw,
			fname.stem = 'raw-cna-calls'
			);
		}

	# rounded CNA calls
	if (! is.null(cna.rounded)) {
		flog.info("Plotting 'rounded CNA calls' density plots..");
		NanoStringNormCNV::make.cna.densities.plots(
			nano.cnas = cna.rounded,
			fname.stem = 'rounded-cna-calls'
			);
		}

	##################################
	### Heatmaps (replicates only) ###
	##################################

	if (! is.null(replicate.eval)) {
		# raw NanoString counts
		flog.info("Plotting raw counts heatmap for replicates..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = replicate.eval$raw.counts,
			fname.stem = 'replicate_raw',
			covs.rows = sample.covs,
			covs.cols = gene.covs.pre.norm,
			);
		
		# normalized NanoString counts
		flog.info("Plotting normalized counts heatmap for replicates..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = replicate.eval$norm.counts,
			fname.stem = 'replicate_norm',
			covs.cols = gene.covs.post.norm,
			covs.rows = sample.covs
			);
		
		# CNA calls
		flog.info("Plotting CNA calls heatmap for replicates..");
		NanoStringNormCNV::make.cna.heatmap(
			nano.cnas = replicate.eval$cna.calls,
			fname.stem = 'replicate_cna-calls',
			rounded = TRUE,
			covs.cols = gene.covs.post.norm,
			covs.rows = sample.covs,
			min.cn = 0
			);

		# replicate CNA call concordance
		flog.info("Plotting replicate CNA call concordance heatmap..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = replicate.eval$concordance,
			fname.stem = 'replicate_cna-concordance',
			covs.cols = gene.covs.post.norm,
			print.ylab = TRUE
			);

		# normalized count correlations
		flog.info("Plotting normalized count correlations heatmap for replicates..");
		NanoStringNormCNV::make.sample.correlations.heatmap(
			nano.counts = log10(replicate.eval$norm.counts + 1),
			fname.stem = 'replicate_norm-count',
			covs = sample.covs
			);


		# normalized count correlations (tumour only)
		flog.info("Plotting normalized count correlations heatmap for replicates (tumour only)..");
		NanoStringNormCNV::make.sample.correlations.heatmap(
			nano.counts = log10(replicate.eval$norm.counts[, which(replicate.eval$count.pheno$Type == 'Tumour')] + 1),
			fname.stem = 'replicate_tumour-only_norm-count',
			covs = sample.covs
			);
		}
	}
