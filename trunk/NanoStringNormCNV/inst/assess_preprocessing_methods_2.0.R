### assess_preprocessing_methods.R #################################################################
# Emilie's analysis modified for new NanoString dataset.

### NOTES ##########################################################################################

### PREAMBLE #######################################################################################
{
	library(mclust);
	library(NanoStringNorm);
	# library(BoutrosLab.utilities.copynumber);
	# library(BoutrosLab.pipeline.limma);
	library(BoutrosLab.plotting.general);
	library(futile.logger);
	library(vsn);
	library(reshape2);
	library(devtools)
	library(getopt);
	# library(BoutrosLab.dist.overload);
	# load_all("~/svn/Resources/code/R/prostate.acgh.biomarkers");
	# source("~/svn/Training/elalonde/OncoScan_reprocess/cna.plotting.functions.R");
	# source("~/svn/Resources/code/R/ParameterEval/R/generate.covariates.R")
	# source("~/svn/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R")
	source("~/svn/Collaborators/RobBristow/nanostring_validation/normalization/call_signature_pga.R")
	load_all("~/svn/Resources/code/R/NanoStringNormCNV/trunk/NanoStringNormCNV");

	# specifically samples with low restriction frag ratios
	dropoutliers <- 0;

	set.seed(12345);

	### FUNCTIONS ######################################################################################
	make_dir_name <- function(params) {
		if (params$perchip == 1) {
			name <- 'perchip';
		} else {
			name <- 'global';
			}

		if (params$ccn == 1) {
			name <- paste(name, 'ccn-sum', sep = '_');
		} else if (params$ccn == 2) {
			name <- paste(name, 'ccn-gm', sep = '_');
			}

		if (params$bc == 1) {
			name <- paste(name, 'bc-mean', sep = '_');
		} else if (params$bc == 2) {
			name <- paste(name, 'bc-m2sd', sep = '_');
		} else if (params$bc == 3) {
			name <- paste(name, 'bc-max', sep = '_');
			}

		if (params$scc == 1) {
			name <- paste(name, 'scc-hk', sep = '_');
		} else if (params$scc == 2) {
			name <- paste(name, 'scc-sum', sep = '_');
		} else if (params$scc == 3) {
			name <- paste(name, 'scc-top-gm', sep = '_');
		} else if (params$scc == 4) {
			name <- paste(name, 'scc-low-gm', sep = '_');
			}

		if (params$inv == 1) {
			name <- paste(name, 'inv', sep = '_');
			}

		if (params$oth == 1) {
			name <- paste(name, 'oth-vsn', sep = '_');
		} else if (params$oth == 2) {
			name <- paste(name, 'oth-quant', sep = '_');
		} else if (params$oth == 3) {
			name <- paste(name, 'oth-rank', sep = '_');
			}

		if (params$matched == 1) {
			name <- paste(name, 'matchedRef', sep = '_');
		} else if (params$matched == 0) {
			name <- paste(name, 'pooledRef', sep = '_');
			}

		if (params$kd == 0) {
			name <- paste(name, 'kd-NS', sep = '_');
		} else if (params$kd == 1) {
			name <- paste(name, 'kd-minMax', sep = '_');
		} else if (params$kd == 2) {
			name <- paste(name, 'kd-default', sep = '_');
		} else if (params$kd == 3) {
			name <- paste(name, 'kd-other', sep = '_');
			}

		if (params$col == 1) {
			name <- paste(name, 'collapsed_by_gene', sep = '_');
			}

		return(name);
		}

	prep_analysis_dir <- function(dir.name, stats = TRUE, plots = TRUE, others = NULL) {
		print(dir("."));
		print(dir.name);
		# make main analysis directory
		if(! file.exists(file.path(dir.name))){
			dir.create(file.path(dir.name), recursive = TRUE);# added recursive arg
			}

		if(stats){
			if(! file.exists(file.path(file.path(dir.name, '/stats/')))){
				dir.create(file.path(file.path(dir.name, '/stats/')));
				}
			}
		if(plots){
			if(! file.exists(file.path(file.path(dir.name, '/plots/')))){
				dir.create(file.path(file.path(dir.name, '/plots/')));
				}
			}

		# if any other directories are requested, make those
		if(! is.null(others)){
			for(n in 1:length(others)){
				if(! file.exists(file.path(file.path(dir.name, '/', others[n])))){
					dir.create(file.path(file.path(dir.name, '/', others[n])));
					}
				}
			}
		}

	check.sample.order <- function(names1, names2) {
		if(all(names1 == names2)){
			return(TRUE);
		}else{
			print(cbind(names1,names2));
			return(FALSE);
			}
		}

	# evaluation.replicates <- function(norm.data, full.phenodata, full.cnas) {
	# 	# prepare data
	# 	nano.reps <- norm.data[, full.phenodata$SampleID[which(full.phenodata$has.repl==1)]];
	# 	pheno.reps <- full.phenodata[which(full.phenodata$has.repl==1), ];
	# 	if(! check.sample.order(colnames(nano.reps), pheno.reps$SampleID)){
	# 		stop("Sorry, samples don't match as expected. See above");
	# 		}

	# 	# order according to samples
	# 	reps.order <- order(pheno.reps$SampleID);
	# 	nano.reps <- nano.reps[, reps.order];
	# 	pheno.reps <- pheno.reps[reps.order , ];
	# 	if(! check.sample.order(colnames(nano.reps), pheno.reps$SampleID)){
	# 		stop("Sorry, samples don't match as expected. See above");
	# 		}
	# 	cna.reps <- full.cnas[ , pheno.reps$SampleID];

	# 	# deal with patient/sample/repl IDs
	# 	samples.backup 		<- pheno.reps$SampleID;
	# 	pheno.reps$replID 	<- pheno.reps$SampleID;
	# 	# pheno.reps$SampleID <- pheno.reps$Patient;
	# 	colnames(cna.reps) 	<- pheno.reps$SampleID;
	# 	pheno.reps$Name 	<- substr(pheno.reps$replID,1,10);
	# 	# names(nano.reps)	<- substr(names(nano.reps),1,10);

	# 	# check variance
	# 	var.matrix <- calculate.replicate.variance(nano.reps, pheno.reps, norm.data$Name);

	# 	### check CNA concordance
	# 	pheno.reps$Name <- substr(pheno.reps$replID,1,10);
	# 	pheno.reps$SampleID <- pheno.reps$Name;
	# 	colnames(cna.reps) <- pheno.reps$Name;
	# 	print(dim(cna.reps));
	# 	print(dim(pheno.reps));
	# 	print(dim(norm.data));
	# 	print(pheno.reps$Name);
	# 	conc.matrix <- calculate.replicate.concordance(cna.reps, pheno.reps, norm.data$Name);
	# 	conc.summary <- apply(conc.matrix, 2, sum)/nrow(conc.matrix);

	# 	pheno.reps$SampleID <- samples.backup;
	# 	return(list(variance = var.matrix, concordance = conc.matrix, conc.summary = conc.summary, pheno = pheno.reps, cnas = cna.reps));
	# 	}

	### SET PARAMETERS #################################################################################
	if (interactive()) {
		opts <- list();
		opts$perchip <- 0;
		opts$ccn  	 <- 2;
		opts$bc 	 <- 2;
		opts$scc 	 <- 2;
		opts$inv 	 <- 1;
		opts$oth 	 <- 0;
		opts$matched <- 0;
		opts$kd 	 <- 3;
		opts$col 	 <- 0;
	} else {
		params <- matrix(
			c(
				'perchip', 'c', 1, 'numeric',
				'ccn', 	   'n', 1, 'numeric',
				'bc', 	   'b', 1, 'numeric',
				'scc', 	   's', 1, 'numeric',
				'inv',	   'i', 1, 'numeric',
				'oth',     'h', 1, 'numeric',
				'matched', 'r', 1, 'numeric',
				'kd',      'k', 1, 'numeric',
				'col',     'o', 1, 'numeric'
				),
			  ncol = 4,
			  byrow = TRUE
	  		);
		opts <- getopt(params);
		}

	# verify arguments
	if(is.null(opts$perchip)) { cat(usage()); q(status = 1) }
	if(is.null(opts$ccn)) 	  { cat(usage()); q(status = 1) }
	if(is.null(opts$bc)) 	  { cat(usage()); q(status = 1) }
	if(is.null(opts$scc)) 	  { cat(usage()); q(status = 1) }
	if(is.null(opts$inv))	  { cat(usage()); q(status = 1) }
	if(is.null(opts$oth)) 	  { cat(usage()); q(status = 1) }
	if(is.null(opts$matched)) { cat(usage()); q(status = 1) }

	home.dir <- make_dir_name(opts);

	# list directories
	raw.dir 	 <- '/.mounts/labs/cpcgene/private/NanoString/nanostring_cancer_panel/nanostring_data/raw_data';
	root.dir     <- '/.mounts/labs/boutroslab/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV';
	training.dir <- '~/svn/Training/Dorota\ Sendorek/NanoStringNormCNV/';

	setwd(root.dir);
	data.dir <- paste0(root.dir, '/test_data/');

	if (dropoutliers == 1) {
		prep_analysis_dir(
			dir.name = paste0('normalization_assessment_outliers_removed/', home.dir),
			stats = FALSE
			);
		out.dir <- paste0(root.dir, '/normalization_assessment_outliers_removed/', home.dir);
	} else {
		prep_analysis_dir(
			dir.name = paste0('normalization_assessment/', home.dir),
			stats = FALSE
			);
		out.dir <- paste0(root.dir, '/normalization_assessment/', home.dir);
		}

	plot.dir <- paste0(out.dir, '/plots/');

	### READ DATA ######################################################################################
	setwd(data.dir);

	# read in raw data, deal with double header (sample name, fragmentation method)
	nano.raw <- read.table('NSrawdata.txt', sep = "\t", skip = 1, stringsAsFactors = FALSE);

	# extract fragmentation methods from raw data (may not be necessary)
	frag.method <- unlist(nano.raw[1,])[-(1:2)];
	nano.raw <- nano.raw[-1,];

	# match Emilie's original NS data format
	names(nano.raw)[1:2] <- c("Accession", "Name");

	nano.raw$CodeClass <- 'Endogenous';
	nano.raw[grep("^POS_[A-Z]", nano.raw$Name),]$CodeClass <- "Positive";
	nano.raw[grep("^NEG_[A-Z]", nano.raw$Name),]$CodeClass <- "Negative";
	nano.raw[grep("^RESTRICTIONSITE", nano.raw$Name),]$CodeClass <- "RestrictionSite";
	nano.raw[grep("_INVCONTROL", nano.raw$Accession),]$CodeClass <- "Invariant";

	# get columns in right order
	nano.raw <- nano.raw[,c(ncol(nano.raw), 2, 1, 3:(ncol(nano.raw) - 1))];

	# read in sample names separately
	sample.names <- scan('NSrawdata.txt', nlines = 1, sep = "\t", what = character());
	sample.names <- sample.names[-(1:2)];

	# set names for consistency
	names(nano.raw)[4:ncol(nano.raw)] <- sample.names;
	names(frag.method) <- sample.names;

	# removing Roche samples for now
	nano.raw <- nano.raw[,-grep("RocheRef", names(nano.raw), perl = TRUE)];

	# modify sample names
	names(nano.raw) <- gsub(" |-",  ".", names(nano.raw));
	names(nano.raw) <- gsub(".RCC", "",  names(nano.raw));

	# integerize it!
	nano.raw[,4:ncol(nano.raw)] <- apply(nano.raw[,4:ncol(nano.raw)], 2, as.integer);

	# fix chr X and Y 'Accession' so they are correctly identified as sex chr and don't get collapsed as single region
	chrXY <- grep("chr[XY]", nano.raw$Accession);
	nano.raw$Name[chrXY] <- paste0("chr", nano.raw$Name[chrXY]);
	nano.raw$Accession[chrXY] <- unlist(lapply(strsplit(nano.raw$Name[chrXY], split = '-'), function(x) x[1]));

	# get phenodata and match Emilie's formatting
	phenodata <- load.phenodata(fname = paste0(data.dir, "/NSannotation.csv"));
	phenodata$Name <- gsub("(.*)\\.M.*", "\\1", phenodata$Name);
	phenodata <- phenodata[, c("SampleID", "Patient", "Name", "cartridge", "type", "ref.name", "has.repl", "sex")];

	# match raw colnames to pheno Sample ID
	check.names <- gsub("_[0-9]+", "", names(nano.raw)[-(1:3)]);
	check.names <- gsub("\\.", "", check.names);
	check.names <- matrix(unlist(strsplit(check.names, "M")), ncol = 2, byrow = T);
	for (i in 1:nrow(check.names)) {
		if (grepl("M[12]", phenodata$SampleID[i])) {
			if (paste0(check.names[i,1], ".M", check.names[i,2]) != phenodata$SampleID[i]) stop("Sample order does not match!");
		} else {
			if (check.names[i,1] != phenodata$SampleID[i]) stop("Sample order does not match!");
			}
		}

	colnames(nano.raw)[-c(1:3)] <- phenodata$SampleID;

	if (! check.sample.order(phenodata$SampleID, colnames(nano.raw)[-c(1:3)])) {
		stop("Sorry, sample order doesn't match after re-naming.");
		}

	# Remove bad normals (low restriction frag ratios) from phenodata
	if (dropoutliers == 1) {
		bad.samples <- read.delim(
			"../normalization_assessment/restriction-fragmentation_low-ratio.txt",
			stringsAsFactors = FALSE,
			header = FALSE
			);
		bad.samples <- bad.samples[,1];

		phenodata <- phenodata[!(phenodata$SampleID %in% bad.samples),];
		nano.raw  <- nano.raw[,!(colnames(nano.raw) %in% bad.samples)];

		# check if any replicates are left
		repls <- phenodata[phenodata$has.repl == 1,];
		repls <- unique(repls[duplicated(repls$Patient) | duplicated(repls$Patient, fromLast = TRUE),]$Patient);
		unlist(lapply(repls, function(x) { any(summary(factor(phenodata[phenodata$Patient == x,]$type)) > 1) }));

		# check that there are no missing references
		if (any(grepl("CPCG", phenodata$ref.name[!(phenodata$ref.name %in% phenodata$SampleID)]))) {
			stop("Sorry, reference samples are missing!");
			}

		# check sample order
		if (! check.sample.order(phenodata$SampleID, colnames(nano.raw)[-c(1:3)])) {
			stop("Sorry, sample order doesn't match after handling outliers.");
			}
		}

	# identifying housekeeping genes (all missing but one)
	nano.raw[nano.raw$Accession %in% qw("ZDHHC5 KIF27 MAGI3 PCDHA9 CPM TMX1 E2F6"), 'CodeClass'] <- 'Housekeeping';

	### Simulating 3 new housekeeping genes by adding noise to original
	original.hk <- nano.raw[nano.raw$CodeClass == 'Housekeeping',];
	simulate.hk <- matrix(nrow = 9, ncol = ncol(nano.raw), dimnames = list(NULL, colnames(nano.raw)));
	simulate.hk <- as.data.frame(simulate.hk);

	simulate.hk$CodeClass <- 'Housekeeping';
	simulate.hk$Accession <- sort(rep(paste0("SIM", 1:3), 3));
	simulate.hk$Name 	  <- as.vector(sapply(paste0("SIM", 1:3, "-"), function(x) paste0(x, 1:3)));

	# simulated HK genes also contain 3 probes
	for (j in 4:ncol(original.hk)) {
		for (i in 1:nrow(original.hk)) {
			# Step 1:
			# adding random noise to original probe at given patient (separately for all 3 simulations)
			noisy.probes <- original.hk[i,j] + rnorm(n = 3, mean = 0, sd = 15);

			# Step 2:
			# shifting mean of probe values per simulation by adding some randomly chosen constant
			noisy.probes <- noisy.probes - c(10, 20, -15);

			# make sure there are no negative values before adding
			noisy.probes[noisy.probes < 1] <- 1;
			simulate.hk[c(i, i + 3, i + 6), j] <- noisy.probes;
			}
		}

	nano.raw <- rbind(nano.raw, simulate.hk);

	# fix gene names to prevent NSN crashing
	nano.raw$Name <- unlist(lapply(strsplit(x = nano.raw$Name, '\\|'), function(f) f[[1]][1]));
	nano.raw$Name <- gsub(x = nano.raw$Name, pattern = '\\.', '');

	# prepare covariates to assess batch effects in NSN
	cartridge.n <- unique(phenodata$cartridge);
	cartridge.matrix <- matrix(nrow = nrow(phenodata), ncol = length(cartridge.n), 1);

	for (n in 1:nrow(phenodata)) {
		cartridge.matrix[n, which(cartridge.n == phenodata$cartridge[n])] <- 2;
		}

	pheno.df <- as.data.frame(cartridge.matrix);
	colnames(pheno.df) <- paste0("Cartridge", seq(1:ncol(cartridge.matrix)));
	pheno.df$Type <- ifelse(phenodata$type == 'Reference', 1, 2);
	rownames(pheno.df) <- phenodata$SampleID;

	# set up params for NSN
	if (opts$ccn == 0) cc.val <- 'none';
	if (opts$ccn == 1) cc.val <- 'sum';
	if (opts$ccn == 2) cc.val <- 'geo.mean';
	# cc.val <- ifelse(opts$ccn == 1, 'geo.mean', 'none');

	if (opts$bc == 0) bc.val <- 'none';
	if (opts$bc == 1) bc.val <- 'mean';
	if (opts$bc == 2) bc.val <- 'mean.2sd';
	if (opts$bc == 3) bc.val <- 'max';
	# bc.val <- ifelse(opts$bc == 1, 'mean.2sd', 'none');

	if (opts$inv == 0) do.rcc.inv.norm <- FALSE;
	if (opts$inv == 1) do.rcc.inv.norm <- TRUE;

	if (opts$scc == 0) sc.val <- 'none';
	if (opts$scc == 1) sc.val <- 'housekeeping.geo.mean';
	if (opts$scc == 2) sc.val <- 'total.sum';
	if (opts$scc == 3) sc.val <- 'top.geo.mean';
	if (opts$scc == 4) sc.val <- 'low.cv.geo.mean';

	if (opts$oth == 0) oth.val <- 'none';
	if (opts$oth == 1) oth.val <- 'vsn';
	if (opts$oth == 2) oth.val <- 'rank.normal';
	if (opts$oth == 3) oth.val <- 'quantile';

	do.nsn.norm <- TRUE;

	# set up kd values if required
	# if (opts$kd > 0) { thresh.method <- 'KD'; }# this variable doesn't actually exist anywhere
	if (opts$kd == 0) kd.vals <- NULL; # 'round'; using NS-provided thresholds
	if (opts$kd == 1) kd.vals <- NULL; # 'round'; using min/max seen in normals
	if (opts$kd == 2) kd.vals <- NULL;						# 'KD'; "pkg defaults"   --ToDo
	if (opts$kd == 3) kd.vals <- c(0.998, 0.79, 0.88, 0.989); # 'KD'; "user-provided"  --ToDo
	# if (opts$kd == 3) kd.vals <- c(0.89, 0.69, 0.65, 0.87); # 'KD'; "user-provided"  --ToDo

	### RUN NORMALIZATION ##############################################################################
	setwd(plot.dir);

	### Positive control normalization + plots
	corrs <- positive.control.norm(nano.raw);
	make.positive.control.plot(
		correlations = corrs,
		covs = phenodata[, c('SampleID', 'type', 'cartridge')]
		);

	### Restriction digestion normalization + plots
	restr.frag.norm.output <- restriction.fragmentation.norm(nano.raw);

	# write bad restr dig samples to file
	write.table(
		restr.frag.norm.output,
		file = paste0(root.dir, "/normalization_assessment/restr-frag-norm_output.txt"),
		quote = FALSE,
		sep = "\t"
		);

	if (dropoutliers == 0) {
		write.table(
			rownames(restr.frag.norm.output[restr.frag.norm.output$ratio < 10,]),
			file = paste0(root.dir, "/normalization_assessment/restriction-fragmentation_low-ratio.txt"),
			quote = FALSE,
			sep = "\t",
			row.names = FALSE,
			col.names = FALSE
			);
		}

	### Invariant probe normalization (this is actually in the main norm fcns)
	inv.probe.norm.output <- invariant.probe.norm(nano.raw, phenodata);
	inv.probe.norm.output <- inv.probe.norm.output[inv.probe.norm.output$CodeClass == 'Invariant',];

	# check that there are no samples with mean invariant probe count < 100
	low.count.samples <- names(which(apply(
		X = nano.raw[nano.raw$CodeClass == "Invariant", -(1:3)],
		MARGIN = 2,
		FUN = mean
		) < 100));

	# drop low counts samples
	if (length(low.count.samples) > 0) {
		print(paste0(
			"Removing low invariant count samples: ",
			paste(low.count.samples, collapse = "  ")
			));

		nano.raw <- nano.raw[, !(colnames(nano.raw) %in% low.count.samples)];
		phenodata <- phenodata[!(phenodata$SampleID %in% low.count.samples),];
		pheno.df <- pheno.df[!(rownames(pheno.df) %in% low.count.samples),];
		}

	# changing ref sample for CPCG248-F1 since original was poor quality
	phenodata[phenodata$SampleID == "CPCG0248F1",]$ref.name <- "CPCG0248B.M1";

	# setting has.repl for replicates to 0
	phenodata[phenodata$SampleID %in% c("CPCG0266B.M2", "CPCG0248B.M1"),]$has.repl <- 0;

	### NanoStringNorm
	setwd(out.dir);

	phenodata$outlier <- 0;
	if (opts$perchip == 1) {
		norm.data <- normalize.per.chip(
			phenodata,
			nano.raw,
			cc.val,
			bc.val,
			sc.val,
			oth.val,
			do.nsn.norm,
			do.rcc.inv.norm,
			pheno.df,
			plot.types = qw('cv mean.sd norm.factors missing RNA.estimates positive.controls')
			);
	} else {
		norm.data <- normalize.global(
			nano.raw,
			cc.val,
			bc.val,
			sc.val,
			oth.val,
			do.nsn.norm,
			do.rcc.inv.norm,
			pheno.df,
			plot.types = qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'),
			pheno = phenodata[,c('SampleID', 'type')]
			);

		# when other normalization is normal.rank we get negative values and need to transform
		if(oth.val == 'rank.normal'){
			norm.data[, -c(1:3)] <- (norm.data[, -c(1:3)] - min(norm.data[, -c(1:3)])) + 0.1;
			}
		}

	if (! check.sample.order(phenodata$SampleID, colnames(norm.data)[-c(1:3)])) {
		stop("Sorry, sample order doesn't match after normalization, see above.");
		}

	### Collapse genes per region if requested #########################################################
	if (opts$col == 1) {

		# save annotations as they will disappear after merge
		norm.annot 	  <- norm.data[, colnames(norm.data) %in% c('Accession', 'CodeClass', 'Name')];
		norm.data 	  <- collapse.genes(norm.data[ , !colnames(norm.data) %in% c('CodeClass', 'Name')]);
		matching.inds <- unlist(lapply(norm.data$Name, function(f) which(norm.annot$Accession == f)[1]));
		norm.data 	  <- cbind(
			Name = norm.data$Name,
			norm.annot[matching.inds, qw("Accession CodeClass")],
			norm.data[, !colnames(norm.data) == 'Name']
			);
		
		# save gene info
		gene.info <- norm.data[,1:3];
		colnames(gene.info)[2] <- 'Symbol';

	} else {

		# save gene info
		gene.info <- norm.data[,1:3];
		colnames(gene.info)[3] <- 'Symbol';

		}

	### Call CNAs ######################################################################################
	if (opts$matched == 1) {
		flog.info('Going to call CNAs with matched normals');

		cna.all <- call.cnas.with.matched.normals(
			normalized.data = norm.data, 
			phenodata = phenodata,
			per.chip = opts$perchip,
			call.method = opts$kd,
			kd.values = kd.vals
			);

		cna.raw <- cna.all$raw;
		cna.rounded <- cna.all$rounded;
	} else {
		flog.info('Going to call CNAs with pooled normals');

		cna.all <- call.cnas.with.pooled.normals(
			normalized.data = norm.data,
			phenodata = phenodata,
			per.chip = opts$perchip,
			call.method = opts$kd,
			kd.values = kd.vals
			);

		cna.rounded <- cna.all$rounded;
		cna.raw <- cna.all$raw;
		cna.normals <- cna.all$normals;
		cna.normals.unadj <- cna.all$normals.unadj;
		}

	has.ref <- which(phenodata$type == 'Tumour');

	# sanity check
	if (! check.sample.order(sub(x = phenodata$SampleID[has.ref], pattern = 'outlier', ''), colnames(cna.rounded))) {
		stop("Sorry, sample order doesn't match after normalization, see above.");
		}

	pheno.cna <- phenodata[has.ref , ];

	### Density plots ##################################################################################
	# normal.for.plot <- norm.data[, phenodata[phenodata$type == "Reference",]$SampleID];
	# normal.for.plot <- as.vector(unlist(normal.for.plot));
	# tumour.for.plot <- norm.data[, phenodata[phenodata$type == "Tumour",]$SampleID];
	# tumour.for.plot <- as.vector(unlist(tumour.for.plot));

	# normal.for.plot <- log10(normal.for.plot + 1);
	# tumour.for.plot <- log10(tumour.for.plot + 1);
	normal.for.plot <- as.vector(na.omit(as.vector(unlist(cna.normals.unadj))));
	tumour.for.plot <- as.vector(na.omit(as.vector(unlist(cna.raw))));
	cnas.for.plot   <- as.vector(na.omit(as.vector(unlist(cna.rounded))));

	density.thresh <- 5;
	normal.for.plot[normal.for.plot > density.thresh] <- density.thresh;
	tumour.for.plot[tumour.for.plot > density.thresh] <- density.thresh;
	cnas.for.plot[	  cnas.for.plot > density.thresh] <- density.thresh;

	normal.density <- density(
		normal.for.plot,
		from = min(c(normal.for.plot, tumour.for.plot)),
		to = max(c(normal.for.plot, tumour.for.plot))
		);
	tumour.density <- density(
		tumour.for.plot,
		from = min(c(normal.for.plot, tumour.for.plot)),
		to = max(c(normal.for.plot, tumour.for.plot))
		);
	cnas.density <- density(
		cnas.for.plot,
		from = min(cnas.for.plot),
		to = max(cnas.for.plot)
		);

	# plotting calls
	plot.colours <- default.colours(3);
	create.densityplot(
		list(
			cnas = cnas.for.plot,
			tumour = tumour.for.plot,
			normal = normal.for.plot
			),
		filename = paste0(plot.dir, "/densityplot_comparison_kd-option-", opts$kd, ".tiff"),
		col = plot.colours,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = list(
						points = list(col = plot.colours, fill = plot.colours, pch = 19),
						text = list(lab = c("CNAs", "tumour", "normal"))
						)
					),
				x = 0.75,
				y = 0.85
				)
			),
		ylimits = c(-0.1, max(c(normal.density$y, tumour.density$y, cnas.density$y)) + .5),
		type = c('l', 'g'),
		xgrid.at = seq(-1, 6, 0.1),
		ygrid.at = seq(0, 20, 0.25)
		);

	# plotting raw tumour and normal
	plot.colours <- default.colours(3)[2:3];
	create.densityplot(
		list(
			tumour = tumour.for.plot,
			normal = normal.for.plot
			),
		filename = paste0(plot.dir, "/densityplot_comparison_raw.tiff"),
		col = plot.colours,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = list(
						points = list(col = plot.colours, fill = plot.colours, pch = 19),
						text = list(lab = c("tumour", "normal"))
						)
					),
				x = 0.75,
				y = 0.85
				)
			),
		ylimits = c(-0.1, max(c(normal.density$y, tumour.density$y)) + .25),
		type = c('l', 'g'),
		xgrid.at = seq(-1, 6, 0.1),
		ygrid.at = seq(0, 20, 0.25)
		);

	# for kd option 3
	density.difference <- normal.density$y - tumour.density$y;
	intersection.point <- normal.density$x[which(diff(density.difference > 0) != 0) + 1];

	# chosen from above/plot
	intersection.point <- c(0.3, 1.67, 2.42, 3.57);

	kd.vals <- list();
	kd.vals[[1]] <- 1 - length(which(tumour.for.plot < intersection.point[1])) / length(tumour.for.plot);
	kd.vals[[2]] <- 1 - length(which(tumour.for.plot < intersection.point[2])) / length(tumour.for.plot);
	kd.vals[[3]] <- length(which(tumour.for.plot < intersection.point[3])) / length(tumour.for.plot);
	kd.vals[[4]] <- length(which(tumour.for.plot < intersection.point[4])) / length(tumour.for.plot);
}

### Evaluate replicates ############################################################################
use.genes <- which(nano.raw$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"));
reps <- evaluate.replicates(
	normalized.data = norm.data[use.genes,],
	phenodata = phenodata,
	cnas = cna.rounded
	);

### OUTPUT #########################################################################################
write.table(
	cbind(norm.data[use.genes, 1:3], cna.raw),
	generate.filename('tmr2ref', 'counts', 'txt'),
	sep = "\t",
	quote = FALSE
	);

write.table(
	cbind(norm.data[use.genes, 1:3], cna.rounded),
	generate.filename('tmr2ref', 'rounded_counts', 'txt'),
	sep = "\t",
	quote = FALSE
	);

write.table(
	norm.data,
	generate.filename('normalized', 'counts', 'txt'),
	sep = "\t",
	quote = FALSE
	);

write.table(
	reps$variance,
	generate.filename('replicate_CNAs', 'variance_matrix', 'txt'),
	sep = "\t",
	quote = FALSE
	);

write.table(
	reps$concordance,
	generate.filename('replicate_CNAs', 'concordance_matrix', 'txt'),
	sep = "\t",
	quote = FALSE
	);

write.table(
	reps$conc.summary,
	generate.filename('replicate_CNAs', 'concordance_summary', 'txt'),
	sep = "\t",
	quote = FALSE
	);

## PLOTS ##################################################################
setwd(plot.dir);

### Custom plots
# covariates
sample.cov <- phenodata[, c('SampleID', 'type', 'cartridge')];
gene.cov <- nano.raw[, c('Name', 'CodeClass')];

# plot input
counts.raw <- nano.raw[,-c(1:3)];
row.names(counts.raw) <- nano.raw$Name;

counts.norm <- norm.data[,-c(1:3)];
row.names(counts.norm) <- norm.data$Name;

# make heatmaps of raw vs normalized data
make.counts.heatmap(
	nano.counts = counts.raw,
	fname.stem = 'raw',
	covs.rows = sample.cov,
	covs.cols = gene.cov
	);
make.counts.heatmap(
	nano.counts = counts.norm,
	fname.stem = 'norm',
	covs.rows = sample.cov,
	covs.cols = gene.cov
	);

make.sample.correlations.heatmap(
	nano.counts = log10(counts.raw + 1),
	fname.stem = 'unnormalized',
	covs = sample.cov
	);
make.sample.correlations.heatmap(
	nano.counts = log10(counts.norm + 1),
	fname.stem = 'normalized',
	covs = sample.cov
	);

# positive control plots
corrs <- positive.control.norm(nano.raw);
make.positive.control.plot(
	correlations = corrs,
	covs = NULL,
	print.x.labels = FALSE
	);
make.positive.control.plot(
	correlations = corrs,
	covs = sample.cov
	);

# and for CNAs
cna.raw[cna.raw > 10] <- 10;
make.cna.heatmap(
	nano.cnas = cna.rounded,
	fname.stem = 'tmr2ref_rounded',
	rounded = TRUE,
	covs.cols = gene.cov,
	covs.rows = sample.cov,
	width = 7
	);
make.cna.heatmap(
	nano.cnas = cna.raw,
	fname.stem = 'tmr2ref_raw',
	rounded = TRUE,
	covs.cols = gene.cov,
	covs.rows = sample.cov,
	width = 7
	);

# make the densities of the CNAs for each call type
make.cna.densities.plots(nano.cnas = cna.rounded, fname.stem = 'tmr2ref_rounded');
make.cna.densities.plots(nano.cnas = cna.raw, fname.stem = 'tmr2ref_raw');

## make heatmaps for reps ################################################
phenodata$replID   <- phenodata$SampleID;
ref.inds 		   <- which(phenodata$SampleID %in% reps$pheno$ref.name);
reps$pheno 		   <- rbind(reps$pheno, phenodata[ref.inds,]);
reps$pheno$Patient <- factor(reps$pheno$Patient);
reps$pheno$type    <- factor(reps$pheno$type);
reps$pheno$outlier <- factor(reps$pheno$outlier, levels = c(0,1));
reps$norm.counts   <- norm.data[ , reps$pheno$SampleID];
reps$raw.counts    <- nano.raw[ , reps$pheno$SampleID];
rownames(reps$raw.counts) <- nano.raw$Name;
#reps$raw.counts   <- nano.raw[ , unlist(lapply(reps$pheno$replID, function(f) which(colnames(nano.raw) == f))) ];

{
	# # make new covariates
	# sample.covs <- generate.covariates(
	# 	x = data.frame(
	# 		Patient = reps$pheno[, "Patient"],
	# 		Type = factor(reps$pheno[, 'type'], levels = c('Blood', 'Tumour')),
	# 		Outlier = reps$pheno$outlier
	# 		),
	# 	colour.list = list(
	# 		Patient = default.colours(nlevels(reps$pheno$Patient)),
	# 		Type = colours()[c(507,532)],
	# 		Outlier = c('white', 'black')
	# 		)
	# 	);

	# sample.covs2 <- generate.covariates(
	# 	x = data.frame(
	# 		Patient = reps$pheno[, "Patient"],
	# 		Outlier = reps$pheno$outlier
	# 		),
	# 	colour.list = list(Patient = default.colours(nlevels(reps$pheno$Patient)), Outlier = c('white', 'black'))
	# 	);

	# sample.legend <- list(
	# 	legend = list(
	# 		colours = colours()[c(507,532)],
	# 		labels = c("Blood", "Tumour")
	# 		),
	#     legend = list(
	# 		colours = default.colours(nlevels(reps$pheno$Patient)),
	# 		labels = levels(reps$pheno$Patient)
	# 		),
	# 	legend = list(
	# 		colours = default.colours(nlevels(nano.raw$CodeClass)),
	# 		labels = levels(nano.raw$CodeClass)
	# 		),
	# 	legend = list(
	# 		colours = c('white', 'black'),
	# 		labels = c('Yes', 'No')
	# 		)
	# 	);
}

make.counts.heatmap(
	reps$norm.counts,
	# reps$norm.counts + 1,
	fname.stem = 'replicate_norm_counts',
	covs.cols = gene.cov,
	covs.rows = sample.cov
	);
make.counts.heatmap(
	reps$raw.counts,
	fname.stem = 'replicate_raw_counts',
	covs.rows = sample.cov,
	covs.cols = gene.cov,
	);

######################################################################
make.cna.heatmap(
	reps$cnas,
	# reps$cnas + 1,
	fname.stem = 'replicate_cnas',
	rounded = TRUE,
	covs.cols = gene.cov,
	covs.rows = sample.cov
	);

make.counts.heatmap(
	reps$concordance,
	fname.stem = 'replicate_cna_concordance',
	covs.cols = gene.cov,
	print.ylab = TRUE
	);

make.sample.correlations.heatmap(
	log10(reps$norm.counts+1),
	fname.stem = 'replicate-normalized',
	covs = sample.cov
	);
make.sample.correlations.heatmap(
	log10(reps$norm.counts[, which(reps$pheno$type == 'Tumour')] + 1),
	fname.stem = 'replicate-tumours-normalized',
	covs = sample.cov
	);

### Save items to compare runs ####################################################################
summary.data <- list();

### save options
summary.data$bc 	 <- opts$bc;
summary.data$cc  	 <- opts$ccn;
summary.data$scc 	 <- opts$scc;
summary.data$other 	 <- opts$oth;
summary.data$matched <- opts$matched;
summary.data$perchip <- opts$perchip;
summary.data$cnas 	 <- opts$kd
summary.data$col 	 <- opts$col;

# sanity check
if(! check.sample.order(reps$pheno$SampleID, colnames(reps$norm.counts))){
	stop("Sorry, sample order doesn't match prior to ARI analysis, see above.");
	}

summary.scores <- score.runs(
	replicates = reps,
	normalized = norm.data,
	cnas = cna.rounded,
	sample.annot = phenodata,
	genes = gene.info,
	normals = cna.normals
	);

### ??
plot.val.rates(summary.scores, fname.stem = 'PCR');

summary.data <- c(summary.data, summary.scores);	# verify this
print(melt(summary.data));

### print to file
setwd(out.dir);
write.table(
	melt(summary.data),
	generate.filename('summary', 'statistics', 'txt'),
	sep = "\t",
	quote = FALSE,
	col.names = FALSE,
	row.names = FALSE
	);

### SESSION_INFO ##################################################################################
save.session.profile(generate.filename('NS_norm_eval', 'Session_Info','txt'), FALSE);

