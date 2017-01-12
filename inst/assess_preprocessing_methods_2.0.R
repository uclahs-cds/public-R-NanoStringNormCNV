### assess_preprocessing_methods_2.0.R #############################################################

### NOTES ##########################################################################################
# Emilie's analysis modified for new NanoString dataset.

### PREAMBLE #######################################################################################
library(mclust);
library(NanoStringNorm);
library(BoutrosLab.plotting.general);
library(futile.logger);
library(reshape2);
library(devtools)
library(getopt);
load_all("~/svn/Resources/code/R/NanoStringNormCNV/trunk/NanoStringNormCNV");

# set options
dropoutliers <- 0;
writetables <- 0;
plotnorm <- 0;

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

	if (params$cnas == 0) {
		name <- paste(name, 'cnas-NS', sep = '_');
	} else if (params$cnas == 1) {
		name <- paste(name, 'cnas-minMax', sep = '_');
	} else if (params$cnas == 2) {
		name <- paste(name, 'cnas-kdDefault', sep = '_');
	} else if (params$cnas == 3) {
		name <- paste(name, 'cnas-kdUser', sep = '_');
		}

	if (params$col == 1) {
		name <- paste(name, 'collapsed_by_gene', sep = '_');
		}

	if (params$vis == 1) {
		name <- paste(name, 'visualized', sep = '_');
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

### SET PARAMETERS #################################################################################
if (interactive()) {
	opts <- list();
	opts$perchip <- 1;
	opts$ccn  	 <- 1;
	opts$bc 	 <- 1;
	opts$scc 	 <- 1;
	opts$inv 	 <- 1;
	opts$oth 	 <- 0;
	opts$matched <- 1;
	opts$cnas 	 <- 1;
	opts$col 	 <- 0;
	opts$vis 	 <- 0;
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
			'cnas',    'k', 1, 'numeric',
			'col',     'o', 1, 'numeric',
			'vis',	   'v', 1, 'numeric'
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
if(is.null(opts$vis))	  { cat(usage()); q(status = 1) }

home.dir <- make_dir_name(opts);
root.dir <- '/.mounts/labs/boutroslab/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV';

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
phenodata$Name[phenodata$Name == "CPCG0346B.M2.01"] <- "CPCG0346.B.M2.01";
phenodata$Name <- gsub("(.*)\\.M.*", "\\1", phenodata$Name);
phenodata <- phenodata[, c("SampleID", "Patient", "Name", "Cartridge", "Type", "ReferenceID", "HasReplicate", "Sex")];

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
	repls <- phenodata[phenodata$HasReplicate == 1,];
	repls <- unique(repls[duplicated(repls$Patient) | duplicated(repls$Patient, fromLast = TRUE),]$Patient);
	unlist(lapply(repls, function(x) { any(summary(factor(phenodata[phenodata$Patient == x,]$Type)) > 1) }));

	# check that there are no missing references
	if (any(grepl("CPCG", phenodata$ReferenceID[!(phenodata$ReferenceID %in% phenodata$SampleID)]))) {
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

# simulated HK genes also contain 3 probes
for (i in 1:nrow(original.hk)) {
	for (j in 1:3) {
		# create a bunch of randomnicity
		randomness <- rnorm(
			n = ncol(nano.raw) - 3,
			mean = rpois(n = ncol(nano.raw) - 3, lambda = 25),
			sd = rgamma(n = ncol(nano.raw) - 3, shape = 50, scale = 2)
			);
		noisy.counts <- original.hk[i, -(1:3)] + floor(randomness);
		noisy.counts[noisy.counts < 1] <- 1;
		nano.raw <- rbind(
			nano.raw,
			c(
				CodeClass = "Housekeeping",
				Name = paste0("SIM", j, "-", i),
				Accession = paste0("SIM", j),
				noisy.counts
				)
			);
		}
	}

# fix gene names to prevent NSN crashing
nano.raw$Name <- unlist(lapply(strsplit(x = nano.raw$Name, '\\|'), function(f) f[[1]][1]));
nano.raw$Name <- gsub(x = nano.raw$Name, pattern = '\\.', '');

# prepare covariates to assess batch effects in NSN
cartridge.n <- unique(phenodata$Cartridge);
cartridge.matrix <- matrix(nrow = nrow(phenodata), ncol = length(cartridge.n), 1);

for (n in 1:nrow(phenodata)) {
	cartridge.matrix[n, which(cartridge.n == phenodata$Cartridge[n])] <- 2;
	}

pheno.df <- as.data.frame(cartridge.matrix);
colnames(pheno.df) <- paste0("Cartridge", seq(1:ncol(cartridge.matrix)));
pheno.df$Type <- ifelse(phenodata$Type == 'Reference', 1, 2);
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
# if (opts$oth == 4) oth.val <- 'zscore';# ignoring because it outputs NAs

do.nsn.norm <- TRUE;

# set up kd values if required
# if (opts$kd > 0) { thresh.method <- 'KD'; }# this variable doesn't actually exist anywhere
if (opts$cnas == 0) kd.vals <- NULL; # 'round'; using NS-provided thresholds
if (opts$cnas == 1) kd.vals <- NULL; # 'round'; using min/max seen in normals
if (opts$cnas == 2) kd.vals <- NULL;						# 'KD'; "pkg defaults"   --ToDo
# if (opts$cnas == 3) kd.vals <- c(0.998, 0.79, 0.88, 0.989); # 'KD'; "user-provided"  --ToDo
if (opts$cnas == 3) kd.vals <- c(0.98, 0.84, 0.92, 0.97); # 'KD'; "user-provided"  --ToDo
# if (opts$cnas == 3) kd.vals <- c(0.89, 0.69, 0.65, 0.87); # 'KD'; "user-provided"  --ToDo

### RUN NORMALIZATION ##############################################################################
setwd(plot.dir);

### Positive control normalization + plots
corrs <- positive.control.norm(nano.raw);
if (plotnorm == 1) {
	make.positive.control.plot(
		correlations = corrs,
		covs = phenodata[, c('SampleID', 'Type', 'Cartridge')]
		);
	}

### Restriction digestion normalization + plots
if (plotnorm == 1) {
	restr.frag.norm.output <- restriction.fragmentation.norm(nano.raw);

	# write bad restr dig samples to file
	write.table(
		restr.frag.norm.output,
		file = paste0(root.dir, "/normalization_assessment/restr-frag-norm_output.txt"),
		quote = FALSE,
		sep = "\t"
		);
	}

### Invariant probe normalization (this is actually in the main norm fcns)
if (plotnorm == 1) {
	inv.probe.norm.output <- invariant.probe.norm(nano.raw, phenodata);
	inv.probe.norm.output <- inv.probe.norm.output[inv.probe.norm.output$CodeClass == 'Invariant',];
	}

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

	nano.raw  <- nano.raw[, !(colnames(nano.raw) %in% low.count.samples)];
	phenodata <- phenodata[ !(phenodata$SampleID %in% low.count.samples),];
	pheno.df  <- pheno.df[  !(rownames(pheno.df) %in% low.count.samples),];
	}

# changing ref sample for CPCG248-F1 since original was poor quality
phenodata[phenodata$SampleID == "CPCG0248F1",]$ReferenceID <- "CPCG0248B.M1";
phenodata[phenodata$SampleID == "CPCG0233F1",]$ReferenceID <- "CPCG0233B.M1";# because 'M2' ref is weird..

# setting HasReplicate for replicates to 0
phenodata[phenodata$SampleID %in% c("CPCG0266B.M2", "CPCG0248B.M1"),]$HasReplicate <- 0;

### NanoStringNorm
setwd(out.dir);

# phenodata$outlier <- 0;
if (opts$perchip == 1) {
	norm.data <- normalize.per.chip(
		pheno = phenodata,
		raw.data = nano.raw,
		cc = cc.val,
		bc = bc.val,
		sc = sc.val,
		oth = oth.val,
		do.nsn = do.nsn.norm,
		do.rcc.inv = do.rcc.inv.norm,
		covs = pheno.df,
		plot.types = qw('cv mean.sd norm.factors missing RNA.estimates positive.controls')
		);
} else {
	norm.data <- normalize.global(
		raw.data = nano.raw,
		cc = cc.val,
		bc = bc.val,
		sc = sc.val,
		oth = oth.val,
		do.nsn = do.nsn.norm,
		do.rcc.inv = do.rcc.inv.norm,
		covs = pheno.df,
		plot.types = qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'),
		pheno = phenodata[,c('SampleID', 'Type')]
		);
	}

if (! check.sample.order(phenodata$SampleID, colnames(norm.data)[-c(1:3)])) {
	stop("Sorry, sample order doesn't match after normalization, see above.");
	}

### Collapse genes per region if requested #########################################################
if (opts$col == 1) {
	norm.data <- collapse.genes(normalized.data = norm.data);
	}

### Call CNAs ######################################################################################
if (opts$matched == 1) {
	flog.info('Going to call CNAs with matched normals');

	cna.all <- call.cnas.with.matched.normals(
		normalized.data = norm.data, 
		phenodata = phenodata,
		per.chip = opts$perchip,
		call.method = opts$cnas,
		kd.values = kd.vals
		);

	cna.raw <- cna.all$raw;
	cna.rounded <- cna.all$rounded;

	has.ref <- which(phenodata$ReferenceID != "missing" & ! is.na(phenodata$ReferenceID));
} else {
	flog.info('Going to call CNAs with pooled normals');

	cna.all <- call.cnas.with.pooled.normals(
		normalized.data = norm.data,
		phenodata = phenodata,
		per.chip = opts$perchip,
		call.method = opts$cnas,
		kd.values = kd.vals
		);

	cna.rounded <- cna.all$rounded;
	cna.raw <- cna.all$raw;
	cna.normals <- cna.all$normals;
	cna.normals.unadj <- cna.all$normals.unadj;

	has.ref <- which(phenodata$Type == 'Tumour');
	}

# sanity check
if (! check.sample.order(sub(x = phenodata$SampleID[has.ref], pattern = 'outlier', ''), colnames(cna.rounded))) {
	stop("Sorry, sample order doesn't match after normalization, see above.");
	}

pheno.cna <- phenodata[has.ref,];


{### Density plots ##################################################################################
	# # normal.for.plot <- norm.data[, phenodata[phenodata$Type == "Reference",]$SampleID];
	# # normal.for.plot <- as.vector(unlist(normal.for.plot));
	# # tumour.for.plot <- norm.data[, phenodata[phenodata$Type == "Tumour",]$SampleID];
	# # tumour.for.plot <- as.vector(unlist(tumour.for.plot));

	# # normal.for.plot <- log10(normal.for.plot + 1);
	# # tumour.for.plot <- log10(tumour.for.plot + 1);

	# # remove sex probes
	# cna.normals.unadj.autosomes <- cna.normals.unadj[!grepl("chr[xy]", tolower(rownames(cna.normals.unadj))),];
	# cna.raw.autosomes 			<- cna.raw[!grepl("chr[xy]", tolower(rownames(cna.raw))),];
	# cna.rounded.autosomes 		<- cna.rounded[!grepl("chr[xy]", tolower(rownames(cna.rounded))),];

	# normal.for.plot <- as.vector(na.omit(as.vector(unlist(cna.normals.unadj.autosomes))));
	# tumour.for.plot <- as.vector(na.omit(as.vector(unlist(cna.raw.autosomes))));
	# cnas.for.plot   <- as.vector(na.omit(as.vector(unlist(cna.rounded.autosomes))));

	# density.thresh <- 5;
	# normal.for.plot[normal.for.plot > density.thresh] <- density.thresh;
	# tumour.for.plot[tumour.for.plot > density.thresh] <- density.thresh;
	# cnas.for.plot[	  cnas.for.plot > density.thresh] <- density.thresh;

	# normal.density <- density(
	# 	normal.for.plot,
	# 	from = min(c(normal.for.plot, tumour.for.plot)),
	# 	to = max(c(normal.for.plot, tumour.for.plot))
	# 	);
	# tumour.density <- density(
	# 	tumour.for.plot,
	# 	from = min(c(normal.for.plot, tumour.for.plot)),
	# 	to = max(c(normal.for.plot, tumour.for.plot))
	# 	);
	# cnas.density <- density(
	# 	cnas.for.plot,
	# 	from = min(cnas.for.plot),
	# 	to = max(cnas.for.plot)
	# 	);

	# # plotting calls
	# plot.colours <- default.colours(3);
	# create.densityplot(
	# 	list(
	# 		cnas = cnas.for.plot,
	# 		tumour = tumour.for.plot,
	# 		normal = normal.for.plot
	# 		),
	# 	filename = paste0(plot.dir, "../../../plots/densityplot_comparison_cnas-option-", opts$cnas, ".tiff"),
	# 	col = plot.colours,
	# 	legend = list(
	# 		inside = list(
	# 			fun = draw.key,
	# 			args = list(
	# 				key = list(
	# 					points = list(col = plot.colours, fill = plot.colours, pch = 19),
	# 					text = list(lab = c("CNAs", "tumour", "normal"))
	# 					)
	# 				),
	# 			x = 0.75,
	# 			y = 0.85
	# 			)
	# 		),
	# 	ylimits = c(-0.1, max(c(normal.density$y, tumour.density$y, cnas.density$y)) + .5),
	# 	type = c('l', 'g'),
	# 	xgrid.at = seq(-1, 6, 0.1),
	# 	ygrid.at = seq(0, 20, 0.25)
	# 	);

	# # plotting raw tumour and normal
	# plot.colours <- default.colours(3)[2:3];
	# create.densityplot(
	# 	list(
	# 		tumour = tumour.for.plot,
	# 		normal = normal.for.plot
	# 		),
	# 	filename = paste0(plot.dir, "../../../plots/densityplot_comparison_raw.tiff"),
	# 	col = plot.colours,
	# 	legend = list(
	# 		inside = list(
	# 			fun = draw.key,
	# 			args = list(
	# 				key = list(
	# 					points = list(col = plot.colours, fill = plot.colours, pch = 19),
	# 					text = list(lab = c("tumour", "normal"))
	# 					)
	# 				),
	# 			x = 0.75,
	# 			y = 0.85
	# 			)
	# 		),
	# 	ylimits = c(-0.1, max(c(normal.density$y, tumour.density$y)) + .25),
	# 	type = c('l', 'g'),
	# 	xgrid.at = seq(-1, 6, 0.1),
	# 	ygrid.at = seq(0, 20, 0.25)
	# 	);

	# # for cnas option 3
	# density.difference <- normal.density$y - tumour.density$y;
	# intersection.point <- normal.density$x[which(diff(density.difference > 0) != 0) + 1];

	# # chosen from above/plot
	# intersection.point <- c(0.3, 1.67, 2.42, 3.57);

	# kd.vals <- list();
	# kd.vals[[1]] <- 1 - length(which(tumour.for.plot < intersection.point[1])) / length(tumour.for.plot);
	# kd.vals[[2]] <- 1 - length(which(tumour.for.plot < intersection.point[2])) / length(tumour.for.plot);
	# kd.vals[[3]] <- length(which(tumour.for.plot < intersection.point[3])) / length(tumour.for.plot);
	# kd.vals[[4]] <- length(which(tumour.for.plot < intersection.point[4])) / length(tumour.for.plot);
	}

### Evaluate replicates ############################################################################
use.genes <- which(norm.data$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"));
reps <- evaluate.replicates(
	normalized.data = norm.data[use.genes,],
	phenodata = phenodata,
	cna.rounded = cna.rounded
	);

### OUTPUT #########################################################################################
# table required for downstream comparison
write.table(
	cbind(norm.data[use.genes, 1:3], cna.rounded),
	generate.filename('tmr2ref', 'rounded_counts', 'txt'),
	sep = "\t",
	quote = FALSE
	);

# optional tables
if (writetables == 1) {
	write.table(
		cbind(norm.data[use.genes, 1:3], cna.raw),
		generate.filename('tmr2ref', 'counts', 'txt'),
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
	}

## PLOTS ##################################################################
if (opts$vis == 1) {
	setwd(plot.dir);

	visualize.results(
		raw.data = nano.raw,
		normalized.data = norm.data,
		phenodata = phenodata,
		cna.rounded = cna.rounded,
		cna.raw = cna.raw,
		replicate.eval = reps,
		max.cn = 10
		);
	}

### Save items to compare runs ####################################################################
summary.data <- list();

### save options
summary.data$bc 	 <- opts$bc;
summary.data$ccn  	 <- opts$ccn;
summary.data$scc 	 <- opts$scc;
summary.data$oth 	 <- opts$oth;
summary.data$matched <- opts$matched;
summary.data$perchip <- opts$perchip;
summary.data$cnas 	 <- opts$cnas;
summary.data$col 	 <- opts$col;
summary.data$inv 	 <- opts$inv;

# sanity check
if(! check.sample.order(reps$count.pheno$SampleID, colnames(reps$norm.counts))){
	stop("Sorry, sample order doesn't match prior to ARI analysis, see above.");
	}

if (opts$matched == 0) {
	score.run.normals <- cna.normals;
} else {
	score.run.normals <- NULL;
	}

summary.scores <- score.runs(
	replicate.eval = reps,
	normalized.data = norm.data,
	cna.rounded = cna.rounded,
	phenodata = phenodata,
	cna.normals = score.run.normals
	);

summary.data <- c(summary.data, summary.scores);

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
