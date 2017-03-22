



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
load_all("~/svn/Resources/code/R/prostate.acgh.biomarkers");
load_all("~/svn/Resources/code/R/NanoStringNormCNV/trunk/NanoStringNormCNV/");
source("~/svn/Training/elalonde/OncoScan_reprocess/cna.plotting.functions.R");
source("~/svn/Resources/code/R/ParameterEval/R/generate.covariates.R")
source("~/svn/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R")
source("~/svn/Collaborators/RobBristow/nanostring_validation/normalization/call_signature_pga.R")
# load_all("~/svn/BoutrosLab/Resources/code/R/prostate.acgh.biomarkers");
# load_all("~/svn/BoutrosLab/Resources/code/R/NanoStringNormCNV/");
# source("~/svn/BoutrosLab/Training/elalonde/OncoScan_reprocess/cna.plotting.functions.R");
# source("~/svn/BoutrosLab/Resources/code/R/ParameterEval/R/generate.covariates.R")
# source("~/svn/BoutrosLab/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R")
# source("~/svn/BoutrosLab/Collaborators/RobBristow/nanostring_validation/normalization/call_signature_pga.R")


dropoutliers  <- 0;#DS
# dropoutliers  <- 1;#EL


### FUNCTIONS #################################################################################
make_dir_name <- function(params){
	if(params$perchip == 1){
		name <- 'perchip';
	}else{
		name <- 'global';
		}
	if(params$ccn ==1){
		name <- paste(name, 'ccn', sep = '_');
		}
	if(params$bc == 1){
		name <- paste(name, 'bc', sep = '_');
		}
	if(params$scc == 1){
		name <- paste(name, 'hk-scc', sep = '_');
		}
	else if(params$scc == 2){
		name <- paste(name, 'top-scc', sep = '_');
		#name <- paste(name, 'inv-hk-scc', sep = '_');
		}
	else if(params$scc == 3){
		name <- paste(name, 'inv-scc', sep = '_');
		}
	else if(params$scc == 4){
		name <- paste(name, 'inv-lm-scc', sep = '_');
		}
	if(params$oth == 1){
		name <- paste(name, 'vsn', sep = '_');
		}
	else if(params$oth == 2){
		name <- paste(name, 'quant', sep = '_');
		}
	else if(params$oth == 3){
		name <- paste(name, 'rank', sep = '_');
		}
	if(params$matched == 1){
		name <- paste(name, 'matchedRef', sep = '_');
		}
	else if(params$matched == 0){
		name <- paste(name, 'pooledRef', sep = '_');
		}
	if(params$kd == 0){
		name <- paste(name, 'threshCNAs', sep = '_');
	}else if(params$kd == 1){
		name <- paste(name, 'kdCNAs1', sep = '_');
	}else if(params$kd == 2){
		name <- paste(name, 'kdCNAs2', sep = '_');
	}else if(params$kd == 3){
		name <- paste(name, 'kdCNAs3', sep = '_');
	}else if(params$kd == 4){
		name <- paste(name, 'kdCNAs4', sep = '_');
		}
	if(params$col == 1){
		name  <- paste(name, 'collapsed_by_sig', sep = '_');
		}

	return(name);
	}

### get parameterization options ###############################################################################
if(interactive()){
	#perchip0 <- ccn1 <- bc1 <- scc3 <- other2 <- matched0 <- kd1 <- col0
	opts <- list();
	opts$perchip <- 0;
	opts$ccn <- 1;
	opts$bc <- 3;
	opts$scc <- 1;
	opts$oth <- 0;
	opts$matched <- 0;
	opts$kd <- 0;
	opts$col <- 0;
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

### If scc = 4 (linear model), we can only very specific paramters
if(opts$scc == 4){
	if(any(c(opts$perchip, opts$ccn, opts$bc, opts$matched) == 1) | opts$oth > 0){
		q(save = 'no');
		}
	}

home.dir <- make_dir_name(opts);
# list directories
root.dir <- '/.mounts/labs/boutroslab/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV/';
setwd(root.dir);
prep_analysis_dir(dir.name = paste0('normalization_assessment/', home.dir), stats = F);
data.dir <- '/.mounts/labs/boutroslab/private/Collaborators/RobBristow/cna_biomarkers/validation/4_Nanostring/data/raw/';
# data.dir <- paste0(root.dir, '/test_data');
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

# drop sample with all 0s
nano.raw <- nano.raw[ , which(colnames(nano.raw) != 'X20140807_Boutros3_36_12')];
nano.raw <- nano.raw[ , which(colnames(nano.raw) != 'X20140910_Boutros23_272_08')];
nano.raw <- nano.raw[ , which(colnames(nano.raw) != 'X20140917_Boutros28_326_02')];
# and 339B1?


### get phenodata
setwd(out.dir);
phenodata <- load.phenodata(rm.outliers = ifelse(dropoutliers == 1, T, F));


nano.raw <- nano.raw[, c(1:3, unlist(lapply(phenodata$Name, function(f) which(colnames(nano.raw) == f))))];
# sanity check
if(! check.sample.order(phenodata$Name, colnames(nano.raw)[-c(1:3)])){
	stop("Sorry, sample order doesn't match after reading and parsing data, see above.");
	}

# rename column headers according to CPCG-ID
cpcg.id <- phenodata$SampleID;
colnames(nano.raw)[-c(1:3)] <- cpcg.id;

# get signature info-- only used if collapsed = TRUE
if(opts$col == 1){
	sig.annot <- read.csv("~/isilon/private/Collaborators/RobBristow/cna_biomarkers/validation/4_Nanostring/array_dev/final/2015-01-14_validation_genes_median_cna_size.txt");
#	                        ~/isilon/private/Collaborators/RobBristow/cna_biomarkers/multivariate/random_forest/2013-04-04_updatedFUtime/2013-04-06_intRisk_glm_collapsed_uncollapsedGeneSignature.txt");
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
}else if(opts$scc==4){
	sc.val <- 'none';
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
	#kd.vals <- c(0.95, 0.92,0.92, 0.95);
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
}else if(opts$scc == 4){
	# normalize with invariant probe, just to get normalized data for normal samples
	norm.data.full <- normalize.global(nano.raw, cc.val, bc.val, 'none', oth.val, do.nsn.norm, TRUE, pheno.df, plot.types=qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'));

	# now read in the pre-processed linear model results
	norm.data <- load.lm.results(collapse = opts$col, outliers.rm = ifelse(dropoutliers == 1, T, F), genes = unique(nano.raw$Accession), patients = phenodata$SampleID);	# need to add unique since when we're using collapsed data, there will be multiple entries for the signature regions
	norm.data <- reshape2::dcast(data=norm.data[,qw("SampleID Gene estimated.log2r")], formula = Gene~SampleID, value.var='estimated.log2r');
	norm.data <- cbind(nano.raw[unlist(lapply(norm.data$Gene, function(f) which(nano.raw$Accession == f)[1])), 1:3], norm.data[, -1]);	# choose only the first hit
	phenodata <- phenodata[unlist(lapply(colnames(norm.data)[-c(1:3)], function(f) which(phenodata$SampleID == f) )) , ];
	if(oth.val == 'rank.normal'){
		norm.data[, -c(1:3)] <- norm.data[, -c(1:3)] + min(norm.data[, -c(1:3)] + 0.1);
		norm.data.full[, -c(1:3)] <- norm.data.full[, -c(1:3)] + min(norm.data.full[, -c(1:3)] + 0.1);
		}
}else{
	norm.data <- normalize.global(raw.data = nano.raw, cc = cc.val, bc = bc.val, sc = sc.val, oth = oth.val, do.nsn = do.nsn.norm, do.rcc.inv = do.rcc.inv.norm, covs = pheno.df, phenodata = pheno.df, plot.types = qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'));
	# norm.data <- normalize.global(nano.raw, cc.val, bc.val, sc.val, oth.val, do.nsn.norm, do.rcc.inv.norm, pheno.df, plot.types=qw('cv mean.sd norm.factors missing RNA.estimates positive.controls'));
	
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

	if(opts$scc == 4){	# we already have the log2ratio so we just need to call CNA states
		cna.raw <- norm.data;
		norm.data.lm <- norm.data;
		is.ref <- grep(x=colnames(norm.data.full), 'B1');
		ref.names <- colnames(norm.data.full)[is.ref];

		# call CNAs
		if(opts$kd <= 1){
			cna.rounded <- norm.data.lm[, -c(1:3)];

			# adjust the T/N so that the median is 0
			norm.data.lm[, -c(1:3)] <- apply(norm.data.lm[, -c(1:3)], 2, function(f) f - median(f));
			# get the t/n signal for normal samples

			if(opts$kd == 0){
				thresh <- c(-1.4, -0.4, 0.4, 1.4);
			}else{
				cna.normals.unadj <- call.copy.number.state(norm.data.full[use.genes,c(1:3, is.ref)], ref.names, per.chip = FALSE, chip.info = pheno.normals, thresh.method = 'none', to.log  =T,  adjust = T)[, -c(1:3)];
				thresh <- c(min(cna.normals.unadj), (min(cna.normals.unadj)+0.65), (max(cna.normals.unadj)-0.7), max(cna.normals.unadj)) - 2;
				}
			cna.normals <- call.copy.number.state(norm.data.full[use.genes,c(1:3, is.ref)], ref.names, per.chip = FALSE, chip.info = pheno.normals, cna.thresh = thresh+2, to.log  =T,  adjust = T)[, -c(1:3)];

			# apply thresholds
			tmp.out <- norm.data.lm[ , -c(1:3)];
			cna.rounded[tmp.out <= thresh[1]] <- 0;
			cna.rounded[tmp.out > thresh[1] & tmp.out <= thresh[2]] <- 1;
			cna.rounded[tmp.out > thresh[2] & tmp.out <= thresh[3]] <- 2;
			cna.rounded[tmp.out > thresh[3] & tmp.out <= thresh[4]] <- 3;
			cna.rounded[tmp.out > thresh[4]] <- 4;

		# Use KD approach
		}else{
			# Call CNAs
			cna.rounded <- apply.kd.cna.thresh(norm.data.lm, kd.thresh = kd.vals);

			# copy over pheno data, just need for chip info
			pheno.normals  <-  phenodata[unlist(lapply(substr(ref.names, 1, 8), function(f) which(phenodata$Patient == f)[1])) , ];
#			# Call CNAs in normals
			cna.normals <- call.copy.number.state(norm.data.full[use.genes,c(1:3, is.ref)], ref.names, per.chip = FALSE, chip.info = pheno.normals, thresh.method = 'KD', kd.vals, adjust = T)[,-c(1:3)];
			}

		# save variables for downstram processing/plotting
		is.tmr <- 1:nrow(phenodata);
		norm.data <- norm.data.lm;
	}else{
		is.tmr <- which(phenodata$type == 'Tumour');
		is.ref <- which(phenodata$type == 'Reference');
		cna.raw <- call.copy.number.state(norm.data[use.genes,], phenodata$SampleID[is.ref], per.chip = opts$perchip, chip.info=phenodata, thresh.method = 'none', adjust = T);
		### make an average ref sample to use to call cnas in normals
		norm.data.tmp <- cbind(norm.data[, c(1:3, (is.ref+3))], avg.ref = apply(X = norm.data[, (is.ref+3)], MARGIN = 1, FUN = mean));
		cna.normals.unadj <- call.copy.number.state(
			norm.data.tmp[use.genes, ], 
			'avg.ref', 
			per.chip = opts$perchip, 
			chip.info=phenodata[is.ref,], 
			thresh.method = 'none', 
			adjust = TRUE
			);
		cna.normals.unadj <- cna.normals.unadj[, -c(1:3)];
		
		if(opts$kd <= 1){
			if(opts$kd == 0){	# predefined
				thresh <- c(0.4, 1.4, 2.4, 3.4);	# DEFAULT
			}else{	# based on normal data
				thresh.offset <- diff(range(cna.normals.unadj) * 0.15);
				print(thresh.offset);
				thresh <- c(min(cna.normals.unadj), min(cna.normals.unadj)+thresh.offset, max(cna.normals.unadj)-thresh.offset, max(cna.normals.unadj));	# try based on % range
				#thresh <- c(min(cna.normals.unadj), min(cna.normals.unadj)+0.65, max(cna.normals.unadj)-0.7, max(cna.normals.unadj));
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
if(opts$scc != 4){
	sample.counts.covs <- make.covariates(phenodata, use.type = TRUE);
}else{
	sample.counts.covs <- samples.cna.covs;
	}
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
if(opts$scc != 4){
	make.counts.heatmap(counts.only, fname.stem='counts', covs.col = gene.covs, covs.row = sample.counts.covs, covs.legend = sample.legend, clust.method = 'euclidean');
	make.samples.correlation.heatmap(log10(nano.raw[, -c(1:3)]+1), fname = generate.filename('unnormalized', 'inter_sample_correlation_heatmap', 'png'), covs = sample.counts.covs, covs.legend = sample.legend)
	}
if(opts$scc == 4){
	make.samples.correlation.heatmap(counts.only, fname = generate.filename('normalized', 'inter_sample_correlation_heatmap', 'png'), covs = sample.counts.covs, covs.legend = sample.legend)
}else{
	make.samples.correlation.heatmap(log10(counts.only+1), fname = generate.filename('normalized', 'inter_sample_correlation_heatmap', 'png'), covs = sample.counts.covs, covs.legend = sample.legend)
}

# and for CNAs
cna.rounded2 <- cna.rounded+1;
make.cna.heatmap(cna.rounded2, fname.stem='tmr2ref_rounded', rounded = TRUE,  covs.col = gene.covs, covs.row = samples.cna.covs, covs.legend = sample.legend);
cna.raw[cna.raw > 10] <- 10;
if(opts$scc != 4){
	cna.raw2 <- cna.raw + 1;
	make.cna.heatmap(cna.raw2, fname.stem='tmr2ref', rounded = TRUE,  covs.col = gene.covs, covs.row = samples.cna.covs, covs.legend = sample.legend);
	}
# candidate genes
plot.ternary.candidate.heatmap(cna.rounded-2, gene.info[use.genes,], generate.filename('candidate_genes', 'heatmap', 'png'), draw.rows = FALSE);


# make the densities of the CNAs for each call type
make.gene.densities.plot(cna.rounded, fname.stem='tmr2ref_rounded', xlab=expression('CNA'));
make.gene.densities.plot(cna.raw, fname.stem='tmr2ref', xlab=expression('CNA'));

### make heatmaps for multi-sample patients (not replicates, but het-study samples)
parafin.samples <- unique(phenodata$Patient[phenodata$tissue == 'FFPE']);
het.samples <- phenodata[phenodata$Patient %in% parafin.samples , ];
het.samples$Patient <- factor(het.samples$Patient, levels = unique(het.samples$Patient));
het.cna <- cna.rounded[ , which(colnames(cna.rounded) %in% het.samples$SampleID)];
colnames(het.cna) <- substr(colnames(het.cna), 1, 8);
het.norm <- counts.only[ , which(colnames(counts.only) %in% het.samples$SampleID)];
# make new covariates
sample.covs <- generate.covariates(
	x = data.frame(
		Patient = het.samples[, "Patient"],
		Type = factor(het.samples[, 'type'], levels = c('Blood', 'Tumoour'))
		),
	colour.list = list(Patient = default.colours(nlevels(het.samples$Patient)), Type= colours()[c(507,532)])
	);
sample.legend <- list(
	legend = list(
		colours = colours()[c(507,532)],
		labels = c("Blood", "Tumour")
		),
    legend = list(
		colours = default.colours(nlevels(het.samples$Patient)),
		labels = levels(het.samples$Patient)
		),
	legend = list(
		colours = default.colours(nlevels(nano.raw$CodeClass)),
		labels = levels(nano.raw$CodeClass)
		)
	);

make.cna.heatmap(het.cna+1, fname.stem='hetsamples_tmr2ref_rounded', rounded = TRUE,  covs.col = gene.covs, covs.row = sample.covs, covs.legend = sample.legend);
plot.ternary.candidate.heatmap(het.cna-2, gene.info[use.genes,], generate.filename('hetsamples_candidate_genes', 'heatmap', 'png'), draw.rows = TRUE);
if(opts$scc == 4){
	make.samples.correlation.heatmap(het.norm, fname = generate.filename('hetsamples', 'normalized_inter_sample_correlation_heatmap', 'png'), covs = sample.covs, covs.legend = sample.legend)
}else{
	make.samples.correlation.heatmap(log10(het.norm+1), fname = generate.filename('hetsamples', 'normalized_inter_sample_correlation_heatmap', 'png'), covs = sample.covs, covs.legend = sample.legend)
}


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

if(opts$scc == 4){
	make.cna.heatmap(reps$norm.counts+1, fname.stem='replicate_norm_counts',  covs.col = gene.covs, covs.row = sample.covs, covs.legend = sample.legend);
}else{
	make.counts.heatmap(reps$norm.counts+1, fname.stem='replicate_norm_counts',  covs.col = gene.covs, covs.row = sample.covs, covs.legend = sample.legend);
	make.counts.heatmap(reps$raw.counts, fname.stem='replicate_raw_counts', covs.row = sample.covs,  covs.col = gene.covs, covs.legend = sample.legend);
	}
make.cna.heatmap(reps$cnas + 1, fname.stem='replicate_cnas', rounded = TRUE,  covs.col = gene.covs, covs.row = sample.covs2, covs.legend = sample.legend);
make.counts.heatmap(reps$concordance, fname.stem='replicate_cna_concordance',  covs.col = gene.covs, print.ylab=T, covs.legend = sample.legend);
plot.ternary.candidate.heatmap(reps$cnas-2, gene.info[use.genes,], generate.filename('replicate_candidate_genes', 'heatmap', 'png'));

if(opts$scc == 4){
	make.samples.correlation.heatmap(reps$norm.counts, fname = generate.filename('replicate_tumours', 'normalized_inter_sample_correlation_heatmap', 'png'), covs = sample.covs, covs.legend = sample.legend)
}else{
	make.samples.correlation.heatmap(log10(reps$norm.counts+1), fname = generate.filename('replicate', 'normalized_inter_sample_correlation_heatmap', 'png'), covs = sample.covs, covs.legend = sample.legend)
	make.samples.correlation.heatmap(log10(reps$norm.counts[, which(reps$pheno$type == 'Tumour')]+1), fname = generate.filename('replicate_tumours', 'normalized_inter_sample_correlation_heatmap', 'png'), covs = sample.covs, covs.legend = sample.legend)
}


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
summary.scores <- score.runs(reps, norm.data, cna.rounded, phenodata, gene.info, ensemble = F, collapsed=opts$col, cna.normals);
plot.val.rates(summary.scores, fname.stem = 'PCR');

summary.data <- c(summary.data, summary.scores);	# verify this
print(melt(summary.data));
### print to file
setwd(out.dir);
write.table(melt(summary.data), generate.filename('summary', 'statistics', 'txt'), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE);


### SESSION_INFO ##################################################################################
save.session.profile(generate.filename('NS_norm_eval', 'Session_Info','txt'), FALSE); # 

