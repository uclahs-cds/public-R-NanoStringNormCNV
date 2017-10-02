### R code from vignette source 'NanoStringNormCNV_Introduction.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=100, signif=3, digits=3)
set.seed(0xdada)

## To create bitmap versions of plots with many dots, circumventing
##   Sweave's fig=TRUE mechanism...
##   (pdfs are too large)
openBitmap = function(nm, rows=1, cols=1) {
  png(paste("NSN-", nm, ".png", sep=""), 
       width=600*cols, height=700*rows, pointsize=14)
  par(mfrow=c(rows, cols), cex=2)
}


###################################################
### code chunk number 2: load.package
###################################################
require("NanoStringNormCNV");


###################################################
### code chunk number 3: eg.load.raw.data
###################################################
require('NanoStringNormCNV');

# load raw count example dataset
data("NanoString.DNA.raw");
str(NanoString.DNA.raw);
print(NanoString.DNA.raw[1:6, 1:7]);


###################################################
### code chunk number 4: eg.load.annotation
###################################################
# load annotation example dataset
data("PhenoData");

# optionally, read in annotation file (same information as above) 
PhenoData <- load.phenodata(
	fname = system.file("extdata", "PhenoData.tsv", package = "NanoStringNormCNV"),
	separator = "tab"
	);

print(head(PhenoData));


###################################################
### code chunk number 5: eg.positive.control.qc (eval = FALSE)
###################################################
## # quality control using positive controls
## r.squared <- positive.control.qc(raw.data = NanoString.DNA.raw);
## 
## # plot R squared values
## make.positive.control.plot(correlations = r.squared, covs = PhenoData);


###################################################
### code chunk number 6: eg.restriction.fragmentation.control.qc.alui (eval = FALSE)
###################################################
## # correctly running QC on AluI-digested samples only
## excl.samples <- PhenoData$SampleID[PhenoData$Fragmentation != "AluI"];
## probe.ratios <- restriction.fragmentation.qc(
## 	raw.data = NanoString.DNA.raw[, ! names(NanoString.DNA.raw) %in% excl.samples]
## 	);


###################################################
### code chunk number 7: eg.restriction.fragmentation.control.qc.all (eval = FALSE)
###################################################
## # running QC on all available samples (\textit{i.e.} AluI-digested and sonicated)
## probe.ratios <- restriction.fragmentation.qc(
## 	raw.data = NanoString.DNA.raw
## 	);


###################################################
### code chunk number 8: eg.invariant.control.qc (eval = FALSE)
###################################################
## # plotting invariant probes
## make.invariant.probe.plot(
## 	inv.probe.counts = NanoString.DNA.raw[NanoString.DNA.raw$CodeClass == 'Invariant', -(1:3)],
## 	tissue.type = PhenoData
## 	);


###################################################
### code chunk number 9: eg.remove.low.qual
###################################################
low.quality 		<- c('CPCG0266B.M1', 'CPCG0248B.M2');
NanoString.DNA.raw  <- NanoString.DNA.raw[, !names(NanoString.DNA.raw) %in% low.quality];
PhenoData   		<- PhenoData[!PhenoData$SampleID %in% low.quality,];


###################################################
### code chunk number 10: eg.update.data
###################################################
# update matched normal and replicate information, as necessary
PhenoData[PhenoData$SampleID == 'CPCG0248F1',]$ReferenceID <- 'missing';
PhenoData[PhenoData$SampleID %in% c('CPCG0266B.M2', 'CPCG0248B.M1'),]$HasReplicate <- 0;


###################################################
### code chunk number 11: eg.write.updated.phenodata (eval = FALSE)
###################################################
## # write updates to file
## write.table(x = PhenoData, file = "PhenoData_updated.csv", sep = ",");


###################################################
### code chunk number 12: eg.norm1
###################################################
# example 1
# perform invariant probe normalization only --cartridges combined
NanoString.DNA.norm <- normalize.global(
	raw.data = NanoString.DNA.raw,
	cc = 'none',
	bc = 'none',
	sc = 'none',
	oth = 'none',
	do.rcc.inv = TRUE,
	covs = NA,
	phenodata = PhenoData
	);


###################################################
### code chunk number 13: eg.norm2
###################################################
# example 2
# perform invariant probe normalization only --cartridges individually
NanoString.DNA.norm <- normalize.per.chip(
	raw.data = NanoString.DNA.raw,
	cc = 'none',
	bc = 'none',
	sc = 'none',
	oth = 'none',
	do.rcc.inv = TRUE,
	covs = NA,
	phenodata = PhenoData
	);


###################################################
### code chunk number 14: eg.norm3
###################################################
# example 3
# include covariates for sample cartridge and sample type 
# covariates must be binary as they are passed directly to NanoStringNorm 'traits'
covs <- as.data.frame(matrix(
	1,
	nrow = nrow(PhenoData),
	ncol = length(unique(PhenoData$Cartridge)),
	dimnames = list(
		PhenoData$SampleID,
		paste0("Cartridge", unique(PhenoData$Cartridge))
		)
	));

for (n in 1:nrow(PhenoData)) {
	covs[n, which(unique(PhenoData$Cartridge) == PhenoData$Cartridge[n])] <- 2;
	}

covs$Type <- ifelse(PhenoData$Type == 'Reference', 1, 2);

NanoString.DNA.norm <- normalize.global(
	raw.data = NanoString.DNA.raw,
	cc = 'none',
	bc = 'none',
	sc = 'none',
	oth = 'none',
	do.rcc.inv = TRUE,
	covs = covs,
	phenodata = PhenoData
	);


###################################################
### code chunk number 15: eg.norm4
###################################################
# same as above but per chip
NanoString.DNA.norm <- normalize.per.chip(
	raw.data = NanoString.DNA.raw,
	cc = 'none',
	bc = 'none',
	sc = 'none',
	oth = 'none',
	do.rcc.inv = TRUE,
	covs = covs,
	phenodata = PhenoData
	);


###################################################
### code chunk number 16: eg.write.normalized (eval = FALSE)
###################################################
## # write normalized counts to file
## write.table(x = NanoString.DNA.norm, file = "normalized_counts.csv", sep = ",");


###################################################
### code chunk number 17: eg.collapse.genes
###################################################
NanoString.DNA.norm.col <- collapse.genes(normalized.data = NanoString.DNA.norm);
print(NanoString.DNA.norm.col[1:6, 1:6]);


###################################################
### code chunk number 18: eg.write.collapsed (eval = FALSE)
###################################################
## # write collapsed data to file
## write.table(x = NanoString.DNA.norm.col, file = "normalized_collapsed_counts.csv", sep = ",");


###################################################
### code chunk number 19: call.cnas.matched.ref
###################################################
# Option 1: call using matched normal reference
cnas <- call.cnas.with.matched.normals(
	normalized.data = NanoString.DNA.norm,
	phenodata = PhenoData,
	per.chip = FALSE,
	call.method = 2,
	kd.values = c(0.99, 0.87, 0.89, 0.96),
	use.sex.info = TRUE
	);


###################################################
### code chunk number 20: call.cnas.pooled.ref
###################################################
# Option 2: call using a pooled normals reference
cnas <- call.cnas.with.pooled.normals(
	normalized.data = NanoString.DNA.norm,
	phenodata = PhenoData,
	per.chip = FALSE,
	call.method = 3,
	use.sex.info = TRUE
	);
# Option 3: call using a pooled normals reference
cnas <- call.cnas.with.pooled.normals(
	normalized.data = NanoString.DNA.norm,
	phenodata = PhenoData,
	per.chip = FALSE,
	call.method = 1,
	use.sex.info = TRUE
	);


###################################################
### code chunk number 21: eg.write.cnas (eval = FALSE)
###################################################
## # write CNAs to file
## write.table(x = cnas$rounded, file = "cnas_rounded.csv", sep = ",");


###################################################
### code chunk number 22: eg.eval.reps
###################################################
# if technical replicates are available
evaluation <- evaluate.replicates(
	phenodata = PhenoData,
	normalized.data = NanoString.DNA.norm,
	cna.rounded = cnas$rounded
	);


###################################################
### code chunk number 23: eg.ari1
###################################################
# how well does the data cluster around the patients from which samples were obtained
patient.ari <- get.ari(
	data.to.cluster = evaluation$cna.calls,
	feature = PhenoData[match(colnames(evaluation$cna.calls), PhenoData$SampleID),]$Patient,
	is.discrete = TRUE
	);


###################################################
### code chunk number 24: eg.ari2
###################################################
# how much does the data cluster around the cartridges on which the samples were processed
# log values, if appropriate
if (all(unlist(NanoString.DNA.norm) >= 0)) {
    count.data <- log10(NanoString.DNA.norm[, -c(1:3)] + 1);
} else {
    count.data <- NanoString.DNA.norm[, -c(1:3)];
    }

cartridge.ari <- get.ari(
    data.to.cluster = count.data,
    feature = PhenoData$Cartridge[match(colnames(NanoString.DNA.norm[, -(1:3)]), PhenoData$SampleID)],
    is.discrete = FALSE
    );


###################################################
### code chunk number 25: eg.vis1a (eval = FALSE)
###################################################
## # plot normalized NanoString counts
## make.counts.heatmap(
## 	nano.counts = NanoString.DNA.norm[, -(1:3)],
## 	fname.stem = 'normalized',
## 	covs.rows = PhenoData[, c('SampleID', 'Type', 'Cartridge')],
## 	covs.cols = NanoString.DNA.raw[, c('Name', 'CodeClass')]
## 	);


###################################################
### code chunk number 26: eg.vis1b (eval = FALSE)
###################################################
## # plot raw NanoString counts
## # make sure raw count data frame has gene names for row names!
## NanoString.DNA.formatted <- NanoString.DNA.raw[, -(1:3)];
## rownames(NanoString.DNA.formatted) <- NanoString.DNA.raw$Name;
## 
## make.counts.heatmap(
## 	nano.counts = NanoString.DNA.formatted,
## 	fname.stem = 'raw',
## 	covs.rows = PhenoData[, c('SampleID', 'Type', 'Cartridge')],
## 	covs.cols = NanoString.DNA.raw[, c('Name', 'CodeClass')]
## 	);


###################################################
### code chunk number 27: eg.vis2a (eval = FALSE)
###################################################
## # plot rounded copy number calls
## make.cna.heatmap(
## 	nano.cnas = cnas$rounded,
## 	fname.stem = 'round',
## 	covs.rows = PhenoData[, c('SampleID', 'Type', 'Cartridge')],
## 	covs.cols = NanoString.DNA.raw[, c('Name', 'CodeClass')],
## 	rounded = TRUE
## 	);


###################################################
### code chunk number 28: eg.vis2b (eval = FALSE)
###################################################
## # plot raw (not rounded) copy number calls
## # first, setting max copy number value at 5
## cnas.raw.max5 <- cnas$raw;
## cnas.raw.max5[cnas.raw.max5 > 5] <- 5;
## 
## make.cna.heatmap(
## 	nano.cnas = cnas.raw.max5,
## 	fname.stem = 'raw',
## 	covs.rows = PhenoData[, c('SampleID', 'Type', 'Cartridge')],
## 	covs.cols = NanoString.DNA.raw[, c('Name', 'CodeClass')],
## 	rounded = FALSE
## 	);


###################################################
### code chunk number 29: eg.vis3 (eval = FALSE)
###################################################
## # plot copy number call density for rounded values
## # two plots: per gene and per sample
## make.cna.densities.plots(
## 	nano.cnas = cnas$rounded
## 	);


###################################################
### code chunk number 30: eg.vis4 (eval = FALSE)
###################################################
## # plot raw NanoString count correlations
## make.sample.correlations.heatmap(
## 	nano.counts = NanoString.DNA.formatted,
## 	covs = PhenoData[, c('SampleID', 'Cartridge', 'Type')]
## 	);


###################################################
### code chunk number 31: eg.vis5 (eval = FALSE)
###################################################
## # alternatively, plot all results using wrapper function
## visualize.results(
##     raw.data = NanoString.DNA.raw,
##     normalized.data = NanoString.DNA.norm,
##     phenodata = PhenoData,
##     cna.rounded = cnas$rounded,
##     cna.raw = cnas$raw,
##     replicate.eval = evaluation,
##     max.cn = 5
##     );


