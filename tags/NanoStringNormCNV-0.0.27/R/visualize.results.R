visualize.results <- function(raw.counts, norm.counts, phenodata = NULL, cna.rounded = NULL, cna.raw = NULL, max.cn = 5, replicate.eval = NULL, exclude.covs = FALSE) {
	# set up covariates
	sample.covs <- gene.covs.pre.norm <- gene.covs.post.norm <- NULL;

	if (!exclude.covs) {
		if (! is.null(phenodata)) {
			sample.covs <- phenodata[, colnames(phenodata) %in% c('SampleID', 'type', 'cartridge'), drop = FALSE];
			if (! any(colnames(sample.covs) %in% 'SampleID')) {
				stop("Must include sample IDs in phenodata!");
				}
			if (ncol(sample.covs) < 2) {
				sample.covs <- NULL;
				flog.warn("Excluding sample covariates due to missing information!");
				}
			}
		gene.covs.pre.norm  <- raw.counts[, c('Name', 'CodeClass')];
		gene.covs.post.norm <- norm.counts[, c('Name', 'CodeClass')];
		}

	# re-format plotting data
	cols.to.remove <- c('CodeClass', 'Name', 'Accession');
	
	raw.counts.reformatted  <- raw.counts[,  !(colnames(raw.counts)  %in% cols.to.remove)];
	norm.counts.reformatted <- norm.counts[, !(colnames(norm.counts) %in% cols.to.remove)];
	
	row.names(raw.counts.reformatted)  <- raw.counts$Name;
	row.names(norm.counts.reformatted) <- norm.counts$Name;

	if (!is.null(cna.raw)) { cna.raw[cna.raw > max.cn] <- max.cn; }

	reps$count.pheno$Patient <- factor(reps$count.pheno$Patient);
	reps$count.pheno$type    <- factor(reps$count.pheno$type);
	if (any(names(reps$count.pheno) == 'outlier')) {
		reps$count.pheno$outlier <- factor(reps$count.pheno$outlier, levels = c(0,1));
		}

	reps$raw.counts <- raw.counts[, reps$count.pheno$SampleID];
	rownames(reps$raw.counts) <- raw.counts$Name;

	##############################
	### Heatmaps (all samples) ###
	##############################

	# raw NanoString counts
	flog.info("Plotting raw counts heatmap..");
	NanoStringNormCNV::make.counts.heatmap(
		nano.counts = raw.counts.reformatted,
		fname.stem = 'raw',
		covs.rows = sample.covs,
		covs.cols = gene.covs.pre.norm
		);

	# normalized NanoString counts
	flog.info("Plotting normalized counts heatmap..");
	NanoStringNormCNV::make.counts.heatmap(
		nano.counts = norm.counts.reformatted,
		fname.stem = 'normalized',
		covs.rows = sample.covs,
		covs.cols = gene.covs.post.norm
		);

	# raw count correlations
	flog.info("Plotting raw count correlations heatmap..");
	NanoStringNormCNV::make.sample.correlations.heatmap(
		nano.counts = log10(raw.counts.reformatted + 1),
		fname.stem = 'raw-count',
		covs = sample.covs
		);

	# normalized count correlations
	flog.info("Plotting normalized count correlations heatmap..");
	NanoStringNormCNV::make.sample.correlations.heatmap(
		nano.counts = log10(norm.counts.reformatted + 1),
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
			nano.counts = reps$raw.counts,
			fname.stem = 'replicate_raw',
			covs.rows = sample.covs,
			covs.cols = gene.covs.pre.norm,
			);
		
		# normalized NanoString counts
		flog.info("Plotting normalized counts heatmap for replicates..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = reps$norm.counts,
			fname.stem = 'replicate_norm',
			covs.cols = gene.covs.post.norm,
			covs.rows = sample.covs
			);
		
		# CNA calls
		flog.info("Plotting CNA calls heatmap for replicates..");
		NanoStringNormCNV::make.cna.heatmap(
			nano.cnas = reps$cna.calls,
			fname.stem = 'replicate_cna-calls',
			rounded = TRUE,
			covs.cols = gene.covs.post.norm,
			covs.rows = sample.covs,
			min.cn = 0
			);

		# replicate CNA call concordance
		flog.info("Plotting replicate CNA call concordance heatmap..");
		NanoStringNormCNV::make.counts.heatmap(
			nano.counts = reps$concordance,
			fname.stem = 'replicate_cna-concordance',
			covs.cols = gene.covs.post.norm,
			print.ylab = TRUE
			);

		# normalized count correlations
		flog.info("Plotting normalized count correlations heatmap for replicates..");
		NanoStringNormCNV::make.sample.correlations.heatmap(
			nano.counts = log10(reps$norm.counts + 1),
			fname.stem = 'replicate_norm-count',
			covs = sample.covs
			);


		# normalized count correlations (tumour only)
		flog.info("Plotting normalized count correlations heatmap for replicates (tumour only)..");
		NanoStringNormCNV::make.sample.correlations.heatmap(
			nano.counts = log10(reps$norm.counts[, which(reps$count.pheno$type == 'Tumour')] + 1),
			fname.stem = 'replicate_tumour-only_norm-count',
			covs = sample.covs
			);
		}
	}
