
### PREAMBLE ######################################################################################
library(mclust);
library(NanoStringNorm);
library(BoutrosLab.utilities.copynumber);
library(BoutrosLab.pipeline.limma);
library(BoutrosLab.plotting.general);
library(futile.logger);
library(vsn);
library(reshape2);
library(devtools)
library(getopt);
library(BoutrosLab.dist.overload);
load_all("~/svn/BoutrosLab/Resources/code/R/prostate.acgh.biomarkers");
load_all("~/svn/BoutrosLab/Resources/code/R/NanoStringNormCNV/");
source("~/svn/BoutrosLab/Training/elalonde/OncoScan_reprocess/cna.plotting.functions.R");
source("~/svn/BoutrosLab/Resources/code/R/ParameterEval/R/generate.covariates.R")
source("~/svn/BoutrosLab/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R");	# *** many functions in here may need to be copied to package ***
source("~/svn/BoutrosLab/Collaborators/RobBristow/nanostring_validation/normalization/call_signature_pga.R")


### FUNCTIONS #################################################################################
### Get user-specified options (these are translated to NSN nomenclature on lines 131-154. Perhaps we should just accept NSN nomenclature instead (*** TO DO ***) 
if(interactive()){	# *** this is useful for development / troubleshooting- let's you manually set the command line parameters so that you can run the script interactively
	opts <- list();
	opts$perchip <- 0;	# 0 means don't process per chip (process all samples together), 1 means process samples per cartridge
	opts$ccn <- 0;	# code count norm option: 0 = none, 1 = geo.mean *** implement other options! ***
	opts$bc <- 0;	# background correction option: 0 = none, 1 = mean.2sd *** implement other options! ***
	opts$scc <- 3; # sample content correction: 1 = housekeeping.geo.mean, 2 = top.geo.mean, 3 = invariant probe normalization (via NSNCNV) *** add an option for none?
	opts$oth <- 0;	# other normalization (not sure this works... to be double checked!)
	opts$matched <- 0; # whether (1) or not (0) to use matched normals for reference. If not, a pooled normal reference is used from all available reference samples.
	opts$kd <- 1;	# how to call CNAs from log2ratio values. 0 = NanoString-defined thresholds from documentation, 1-4 = different kernal density (KD) values. 1 = define thresholds based on log2ratio values observed in reference samples (see lines 244-247) and 2-4 are values I defined based on my dataset, and probably not at all generic. *** In reality we should provide 1 set of default KD values, and have another option where users can define them (instead of 3 predefined value sets (kd = 2-4) like we have now). See lines 156-166 for the values currently being used.
	opts$col <- 0;	# whether (1) or not(0) to collapse probes to larger unit (e.g. usually 3 probes per gene)
}else{
	params <- matrix(
		c(
			'perchip', 'c', 1, 'numeric',
			'ccn', 'n', 1, 'numeric',
			'bc', 'b', 1, 'numeric',
			'scc', 's', 1, 'numeric',
			'oth', 'h', 1, 'numeric',
			'matched', 'r', 1, 'numeric',
			'kd', 'k', 1, 'numeric',
			'col', 'o', 1, 'numeric'
			),
		  ncol = 4,
		  byrow = TRUE
  		);
	opts <- getopt(params);
	}
# verify arguments
if(is.null(opts$perchip)) { cat(usage()); q(status = 1) }
if(is.null(opts$ccn)) { cat(usage()); q(status = 1) }
if(is.null(opts$bc)) { cat(usage()); q(status = 1) }
if(is.null(opts$scc)) { cat(usage()); q(status = 1) }
if(is.null(opts$oth)) { cat(usage()); q(status = 1) }
if(is.null(opts$matched)) { cat(usage()); q(status = 1) }


home.dir <- make_dir_name(opts);
# list directories
root.dir <- '/isilon/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV/';
setwd(root.dir);
prep_analysis_dir(dir.name = paste0('normalization_assessment/', home.dir), stats = F);
data.dir <- paste0(root.dir, '/test_data');
out.dir  <- paste0(root.dir, '/normalization_assessment/', home.dir);
plot.dir <- paste0(out.dir, '/plots');

### READ DATA #####################################################################################
# set working directory
setwd(data.dir);

# read in the RCC files
data.raw <- read.markup.RCC(
	rcc.path = data.dir,
	rcc.pattern = "*.RCC"
	);

# only keep the counts and not the header
nano.raw <- data.raw$x;

# drop sample with all 0s *** need to check if any such samples and remove programatically! ***
nano.raw <- nano.raw[ , which(colnames(nano.raw) != 'X20140807_Boutros3_36_12')];
nano.raw <- nano.raw[ , which(colnames(nano.raw) != 'X20140910_Boutros23_272_08')];
nano.raw <- nano.raw[ , which(colnames(nano.raw) != 'X20140917_Boutros28_326_02')];


### get phenodata (*** sample annotations-- will need to specify format of this data in docs so that this function can be published ***)
setwd(out.dir);
phenodata <- load.phenodata(rm.outliers = ifelse(dropoutliers == 1, T, F));


nano.raw <- nano.raw[, c(1:3, unlist(lapply(phenodata$Name, function(f) which(colnames(nano.raw) == f))))];
# sanity check
if(! check.sample.order(phenodata$Name, colnames(nano.raw)[-c(1:3)])){
	stop("Sorry, sample order doesn't match after reading and parsing data, see above.");
	}


# get signature info-- only used if collapsed = TRUE
if(opts$col == 1){
	sig.annot <- read.csv("~/isilon/private/Collaborators/RobBristow/cna_biomarkers/validation/4_Nanostring/array_dev/final/2015-01-14_validation_genes_median_cna_size.txt");
	sig.annot <- merge(sig.annot[sig.annot$assay == 'signature' , ], nano.raw[, 1:3], by.x = 'ns.symbol', by.y = 'Accession', all.x = F, all.y = F);
	in.sig <- unlist(lapply(sig.annot$ns.symbol, function(f) which(nano.raw$Accession == f)));
	nano.raw$Accession[in.sig] <- paste0('sig', sig.annot$sig.row);
	}

# sanity check
if(! check.sample.order(phenodata$SampleID, colnames(nano.raw)[-c(1:3)])){
	stop("Sorry, sample order doesn't match after handling outliers, see above.");
	}

# specify housekeeping genes
nano.raw[nano.raw$Accession %in% qw("ZDHHC5 KIF27 MAGI3 PCDHA9 CPM TMX1 E2F6"), 'CodeClass'] <- 'Housekeeping';

# fix gene names to prevent NSN crashing
nano.raw$Name <- unlist(lapply(strsplit(x=nano.raw$Name, '\\|'), function(f) f[[1]][1]));
nano.raw$Name <- gsub(x=nano.raw$Name, pattern = '\\.', '');

# prepare covariates to assess batch effects in NSN
cartridge.n <- unique(phenodata$cartridge);
cartridge.matrix <- matrix(nrow = nrow(phenodata), ncol = length(cartridge.n), 1);
for(n in 1:nrow(phenodata)){
	cartridge.matrix[n, which(cartridge.n == phenodata$cartridge[n])] <- 2;
	}
pheno.df <- as.data.frame(cartridge.matrix);
colnames(pheno.df) <- paste0("Cartridge", seq(1:ncol(cartridge.matrix)));
pheno.df$Type <- ifelse(phenodata$type == 'Reference', 1, 2);
rownames(pheno.df) <- phenodata$SampleID;

### set up params for NSN
cc.val <- ifelse(opts$ccn ==1, 'geo.mean', 'none');
bc.val <- ifelse(opts$bc ==1, 'mean.2sd', 'none');
do.rcc.inv.norm <- FALSE;
if(opts$scc==1){
	sc.val <- 'housekeeping.geo.mean';
}else if(opts$scc==2){
	sc.val <- 'top.geo.mean';
}else if(opts$scc==3){
	sc.val <- 'none';
	do.rcc.inv.norm <- TRUE;
	}
oth.val <- 'none';
if(opts$oth == 1){
	oth.val <- 'vsn';
}else if(opts$oth == 2){
	oth.val <- 'rank.normal';
}else if(opts$oth == 3){
	oth.val <- 'quantile';
	}
do.nsn.norm <- FALSE;
if(opts$scc < 3 | opts$ccn == 1 | opts$bc == 1 | opts$oth > 0){
	do.nsn.norm <- TRUE;
	}

# set up kd values if required
if(opts$kd > 0){ thresh.method <- 'KD'; }
if(opts$kd == 1){ # try min/max seen in normals
	kd.vals <- c(0.85,0.95);
}else if(opts$kd == 2){
	kd.vals <- c(0.9, 0.87,0.93, 0.96);
}else if(opts$kd == 3){
	kd.vals <- c(0.9, 0.8, 0.87, 0.9);
}else if(opts$kd == 4){
	kd.vals <- c(0.9, 0.885, 0.92, 0.97);
	}

### RUN NORMALIZATION ##############################################################################
phenodata$outlier <- 0;
if(opts$perchip == 1){	# each chip separately
	norm.data <- normalize.per.chip(phenodata, nano.raw, cc.val, bc.val, sc.val, oth.val, do.nsn.norm, do.rcc.inv.norm, pheno.df, plot.types=qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'));
}else{
	norm.data <- normalize.global(nano.raw, cc.val, bc.val, sc.val, oth.val, do.nsn.norm, do.rcc.inv.norm, pheno.df, plot.types=qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'));
	### When other normalization is normal.rank we get negative values and need to transform
	if(oth.val == 'rank.normal'){
		norm.data[, -c(1:3)] <- (norm.data[, -c(1:3)] - min(norm.data[, -c(1:3)])) + 0.1;
		}
	}

# sanity check
if(! check.sample.order(phenodata$SampleID, colnames(norm.data)[-c(1:3)])){
	stop("Sorry, sample order doesn't match after normalization, see above.");
	}

### Collapse genes per region if requested ############################################################
if(opts$col == 1){
	# save annotations as they will dissapear after merge
	norm.annot <- norm.data[, colnames(norm.data) %in% c('Accession', 'CodeClass', 'Name')];
	norm.data <- collapse.genes(norm.data[ , !colnames(norm.data) %in% c('CodeClass', 'Name')]);
	matching.inds <- unlist(lapply(norm.data$Name, function(f) which(norm.annot$Accession == f)[1]));
	norm.data <- cbind(Name = norm.data$Name, norm.annot[matching.inds, qw("Accession CodeClass")], norm.data[, !colnames(norm.data) == 'Name']);
	# save gene info
	gene.info <- norm.data[,1:3];
	colnames(gene.info)[2] <- 'Symbol';
}else{
	# save gene info
	gene.info <- norm.data[,1:3];
	colnames(gene.info)[3] <- 'Symbol';
	}


### Call CNAs #########################################################################################
use.genes <- which(norm.data$CodeClass %in% qw("Endogenous Housekeeping Invariant"));
cna.normals <- matrix(nrow = length(use.genes), ncol = length(which(phenodata$type == 'Reference')));
cna.normals.unadj <- matrix(nrow = length(use.genes), ncol = length(which(phenodata$type == 'Reference')));
if(opts$matched == 1){
	flog.info('Going to call CNAs with matched normals');
	
	#### *** put code below in it's own function! ***
	has.ref <- which(phenodata$ref.name != 'missing' & phenodata$type == 'Tumour');
	cna.raw <- matrix(nrow = length(use.genes), ncol = length(has.ref));
	cna.rounded <- matrix(nrow = length(use.genes), ncol = length(has.ref));

	# iterate through each sample here
	for(tmr in 1:length(has.ref)){
		# find indices for tmr and ref
		tmr.ind <- which(colnames(norm.data) == phenodata$SampleID[has.ref[tmr]]);
		ref.ind <- which(colnames(norm.data) == phenodata$ref.name[has.ref[tmr]]);
		flog.info('The indices for tmr and ref for sample %s: %s and %s', capture = T, c(phenodata$SampleID[has.ref[tmr]], tmr.ind, ref.ind));
		cna.raw[, tmr] <- call.copy.number.state(norm.data[use.genes, c(1:3, tmr.ind, ref.ind)], phenodata$ref.name[has.ref[tmr]], thresh.method = 'none', multi.factor = 2)[,4];
		if(opts$kd <= 1){
			cna.rounded[,tmr] <- call.copy.number.state(norm.data[use.genes, c(1:3, tmr.ind, ref.ind)], phenodata$ref.name[has.ref[tmr]], per.chip = opts$perchip, chip.info = phenodata)[,4];
		}else{
			cna.rounded[,tmr] <- call.copy.number.state(norm.data[use.genes, c(1:3, tmr.ind, ref.ind)], phenodata$ref.name[has.ref[tmr]], per.chip = opts$perchip, chip.info = phenodata, thresh.method = 'KD')[,4];
			cna.normals <- cna.rounded[ , phenodata$SampleID[is.ref]];
			}
		}
}else{ #pooled normals
	flog.info('Going to call CNAs with pooled normals');

	#### *** put code below in it's own function! ***
	
	is.tmr <- which(phenodata$type == 'Tumour');
	is.ref <- which(phenodata$type == 'Reference');
	cna.raw <- call.copy.number.state(norm.data[use.genes,], phenodata$SampleID[is.ref], per.chip = opts$perchip, chip.info=phenodata, thresh.method = 'none', adjust = T);
	### make an average ref sample to use to call cnas in normals
	norm.data.tmp <- cbind(norm.data[, c(1:3, (is.ref+3))], avg.ref = apply(X = norm.data[, (is.ref+3)], MARGIN = 1, FUN = mean));
	cna.normals.unadj <- call.copy.number.state(norm.data.tmp[use.genes, ], 'avg.ref', per.chip = opts$perchip, chip.info=phenodata[is.ref,], thresh.method = 'none', adjust = T);
	cna.normals.unadj <- cna.normals.unadj[, -c(1:3)];
	
	if(opts$kd <= 1){
		if(opts$kd == 0){	# predefined
			thresh <- c(0.4, 1.4, 2.4, 3.4);	# DEFAULT
		}else{	# based on normal data
			thresh.offset <- diff(range(cna.normals.unadj) * 0.15);
			print(thresh.offset);
			thresh <- c(min(cna.normals.unadj), min(cna.normals.unadj)+thresh.offset, max(cna.normals.unadj)-thresh.offset, max(cna.normals.unadj));	# try based on % range
			}
		# call  CNAs for tumours based on derivd thresh (above)
		cna.rounded <- call.copy.number.state(norm.data[use.genes,], phenodata$SampleID[is.ref], per.chip = opts$perchip, chip.info = phenodata, adjust = T, cna.thresh = thresh);
	}else{
		# Call CNAs with KD
		cna.rounded <- call.copy.number.state(norm.data[use.genes,], phenodata$SampleID[is.ref], per.chip = opts$perchip, chip.info = phenodata, thresh.method = 'KD', to.log = FALSE, kd.vals = kd.vals, adjust = T);
		}
	cna.normals <- cna.rounded[ , phenodata$SampleID[is.ref]];
	}


# check for columns with all NA-- need to drop those (happens when perchip is true and matched is false, and there are no reference samples on that chip)
if(any(apply(cna.raw, 2, function(f) all(is.na(f))))){
	all.na <- which(as.vector(apply(cna.raw[,-c(1:3)], 2, function(f) all(is.na(f)))));
	cna.raw <- cna.raw[ , -(all.na+3)];
		}

	cna.raw <- as.matrix(cna.raw[, grep(x=colnames(cna.raw), 'F|P[0-9]')]);
	cna.rounded <- as.matrix(cna.rounded[, grep(x = colnames(cna.rounded), 'F|P[0-9]')]);
	has.ref <- is.tmr[colnames(norm.data)[is.tmr+3] %in% colnames(cna.rounded)];
	}

# sanity check
if(! check.sample.order(sub(x=phenodata$SampleID[has.ref], pattern = 'outlier', ''), colnames(cna.rounded))){
	stop("Sorry, sample order doesn't match after normalization, see above.");
	}
pheno.cna <- phenodata[has.ref , ];
colnames(cna.rounded) <- phenodata$SampleID[has.ref];


### Evaluate replicates ######################################################################
reps <- evaluation.replicates(norm.data[use.genes,], pheno.cna, cna.rounded);

### OUTPUT #################################################################
write.table(cbind(norm.data[use.genes, 1:3], cna.raw), generate.filename('tmr2ref', 'counts', 'txt'), sep = "\t", quote = FALSE);
write.table(cbind(norm.data[use.genes, 1:3], cna.rounded), generate.filename('tmr2ref', 'rounded_counts', 'txt'), sep = "\t", quote = FALSE);
write.table(norm.data, generate.filename('normalized', 'counts', 'txt'), sep = "\t", quote = FALSE);
write.table(reps$concordance, generate.filename('replicate_CNAs', 'concordance_matrix', 'txt'), sep = "\t", quote = FALSE);
write.table(reps$conc.summary, generate.filename('replicate_CNAs', 'concordance_summary', 'txt'), sep = "\t", quote = FALSE);



## PLOTS ##################################################################
setwd(plot.dir);


### Set up covariates and legend for custom plots
nano.raw$CodeClass <- as.factor(nano.raw$CodeClass);
phenodata$type <- as.factor(phenodata$type);
phenodata$Patient <- as.factor(phenodata$Patient);
phenodata$outlier <- as.factor(phenodata$outlier);
pheno.cna$cartridge <- as.factor(pheno.cna$cartridge);
phenodata$cartridge <- as.factor(phenodata$cartridge);
pheno.cna$outlier <- as.factor(pheno.cna$outlier);
samples.cna.covs <- make.covariates(pheno.cna, use.type = F);
sample.counts.covs <- make.covariates(phenodata, use.type = TRUE);
gene.covs <- make.gene.covariates(nano.raw$CodeClass);
sample.legend <- list(
	legend = list(
		colours = colours()[c(507,532)],
		labels = c("Blood", "Tumour")
		),
	legend = list(
		colours = default.colours(nlevels(nano.raw$CodeClass)),
		labels = levels(nano.raw$CodeClass)
		)
	);


### Custom plots
# make heatmaps of raw vs normalized data
counts.only <- norm.data[,-c(1:3)];
make.counts.heatmap(nano.raw[, -c(1:3)], fname.stem='raw',  covs.col = gene.covs, covs.row = sample.counts.covs, covs.legend = sample.legend);
make.counts.heatmap(counts.only, fname.stem='counts', covs.col = gene.covs, covs.row = sample.counts.covs, covs.legend = sample.legend, clust.method = 'euclidean');
make.samples.correlation.heatmap(log10(nano.raw[, -c(1:3)]+1), fname = generate.filename('unnormalized', 'inter_sample_correlation_heatmap', 'png'), covs = sample.counts.covs, covs.legend = sample.legend)
make.samples.correlation.heatmap(log10(counts.only+1), fname = generate.filename('normalized', 'inter_sample_correlation_heatmap', 'png'), covs = sample.counts.covs, covs.legend = sample.legend)

# and for CNAs
cna.rounded2 <- cna.rounded+1;
make.cna.heatmap(cna.rounded2, fname.stem='tmr2ref_rounded', rounded = TRUE,  covs.col = gene.covs, covs.row = samples.cna.covs, covs.legend = sample.legend);
cna.raw[cna.raw > 10] <- 10;
cna.raw2 <- cna.raw + 1;
make.cna.heatmap(cna.raw2, fname.stem='tmr2ref', rounded = TRUE,  covs.col = gene.covs, covs.row = samples.cna.covs, covs.legend = sample.legend);

# candidate genes
plot.ternary.candidate.heatmap(cna.rounded-2, gene.info[use.genes,], generate.filename('candidate_genes', 'heatmap', 'png'), draw.rows = FALSE);


# make the densities of the CNAs for each call type
make.gene.densities.plot(cna.rounded, fname.stem='tmr2ref_rounded', xlab=expression('CNA'));
make.gene.densities.plot(cna.raw, fname.stem='tmr2ref', xlab=expression('CNA'));



### make heatmaps for reps ################################################
phenodata$replID <- phenodata$SampleID;
ref.inds <- which(phenodata$SampleID %in% reps$pheno$ref.name);
reps$pheno <- rbind(reps$pheno, phenodata[ref.inds,]);
reps$pheno$Patient <- factor(reps$pheno$Patient);
reps$pheno$type <- factor(reps$pheno$type);
reps$pheno$outlier <- factor(reps$pheno$outlier, levels = c(0,1));
reps$norm.counts <- norm.data[ , reps$pheno$SampleID];
reps$raw.counts <- nano.raw[ , reps$pheno$SampleID];
#reps$raw.counts <- nano.raw[ , unlist(lapply(reps$pheno$replID, function(f) which(colnames(nano.raw) == f))) ];

# make new covariates
sample.covs <- generate.covariates(
	x = data.frame(
		Patient = reps$pheno[, "Patient"],
		Type = factor(reps$pheno[, 'type'], levels = c('Blood', 'Tumour')),
		Outlier = reps$pheno$outlier
		),
	colour.list = list(Patient = default.colours(nlevels(reps$pheno$Patient)), Type= colours()[c(507,532)], Outlier = c('white', 'black'))
	);
sample.covs2 <- generate.covariates(
	x = data.frame(
		Patient = reps$pheno[, "Patient"],
		Outlier = reps$pheno$outlier
		),
	colour.list = list(Patient = default.colours(nlevels(reps$pheno$Patient)), Outlier = c('white', 'black'))
	);
sample.legend <- list(
	legend = list(
		colours = colours()[c(507,532)],
		labels = c("Blood", "Tumour")
		),
    legend = list(
		colours = default.colours(nlevels(reps$pheno$Patient)),
		labels = levels(reps$pheno$Patient)
		),
	legend = list(
		colours = default.colours(nlevels(nano.raw$CodeClass)),
		labels = levels(nano.raw$CodeClass)
		),
	legend = list(
		colours = c('white', 'black'),
		labels = c('Yes', 'No')
		)
	);

make.counts.heatmap(reps$norm.counts+1, fname.stem='replicate_norm_counts',  covs.col = gene.covs, covs.row = sample.covs, covs.legend = sample.legend);
make.counts.heatmap(reps$raw.counts, fname.stem='replicate_raw_counts', covs.row = sample.covs,  covs.col = gene.covs, covs.legend = sample.legend);
make.cna.heatmap(reps$cnas + 1, fname.stem='replicate_cnas', rounded = TRUE,  covs.col = gene.covs, covs.row = sample.covs2, covs.legend = sample.legend);
make.counts.heatmap(reps$concordance, fname.stem='replicate_cna_concordance',  covs.col = gene.covs, print.ylab=T, covs.legend = sample.legend);
plot.ternary.candidate.heatmap(reps$cnas-2, gene.info[use.genes,], generate.filename('replicate_candidate_genes', 'heatmap', 'png'));
make.samples.correlation.heatmap(reps$norm.counts, fname = generate.filename('replicate_tumours', 'normalized_inter_sample_correlation_heatmap', 'png'), covs = sample.covs, covs.legend = sample.legend)


### Save items to compare runs ####################################################################
summary.data <- list();
### save options
summary.data$bc <- opts$bc;
summary.data$cc <- opts$ccn;
summary.data$scc <- opts$scc;
summary.data$other <- opts$oth;
summary.data$matched <- opts$matched;
summary.data$perchip <- opts$perchip;
summary.data$cnas <- opts$kd
summary.data$col <- opts$col;

# sanity check
if(! check.sample.order(reps$pheno$SampleID, colnames(reps$norm.counts))){
	stop("Sorry, sample order doesn't match prior to ARI analysis, see above.");
	}
if(opts$scc == 4){
	gene.info$Symbol[gene.info$Symbol == 'sig7'] <- 'TERT';
#	summary.scores <- score.runs(reps, norm.data, cna.rounded, phenodata, gene.info, ensemble = T);
}


### *** see score.runs function for the different ways to evaluate the preprocessing method - not all will be applicable to new dataset!
summary.scores <- score.runs(reps, norm.data, cna.rounded, phenodata, gene.info, ensemble = F, collapsed=opts$col, cna.normals);
summary.data <- c(summary.data, summary.scores);
print(melt(summary.data));
### print to file
setwd(out.dir);
write.table(melt(summary.data), generate.filename('summary', 'statistics', 'txt'), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE);


### SESSION_INFO ##################################################################################
save.session.profile(generate.filename('NS_norm_eval', 'Session_Info','txt'), FALSE); # 

