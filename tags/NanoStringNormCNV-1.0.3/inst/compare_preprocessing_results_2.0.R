### compare_preprocessing_results_2.0.R ###########################################################

### PREAMBLE ######################################################################################
library(mclust);
library(NanoStringNorm);
library(BoutrosLab.utilities);
library(BoutrosLab.pipeline.limma);
library(BoutrosLab.plotting.general);
library(futile.logger);
library(reshape2);
library(devtools);
library(getopt);
library(randomForest);
load_all("~/svn/Resources/code/R/prostate.acgh.biomarkers");
load_all("~/svn/Resources/code/R/NanoStringNormCNV/trunk/NanoStringNormCNV");
source("~/svn/Training/elalonde/OncoScan_reprocess/cna.plotting.functions.R");
source("~/svn/Resources/code/R/ParameterEval/R/generate.covariates.R")
source("~/svn/Resources/code/R/BoutrosLab.statistics.general/R/get.pve.R")
source("~/svn/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R")

# set project!
if (interactive()) {
	opts <- list();
	opts$proj.stem <- 'bristow';
#	opts$proj.stem <- 'nsncnv';
	# opts$proj.stem <- 'nsncnv_col';
} else {
	params <- matrix(
		c(
			'proj.stem', 'p', 1, 'character'
			),
		  ncol = 4,
		  byrow = TRUE
  		);
	opts <- getopt(params);
	}
proj.stem <- opts$proj.stem;

remove.any.na.runs <- TRUE;

parameters <- qw('perchip bc ccn scc matched oth cnas col');
score.and.sort <- rbind(
	c('ari.chip', FALSE, 0),
	c('ari.pts.normcor', TRUE, 0),
	c('ari.type', TRUE, 0),
	c('ari.pts', TRUE, 0),
	c('sd.inv', FALSE, 0),
	c('sd.hk', FALSE, 0),
	c('sd.inv.and.hk', FALSE, 0),
	c('replicates.conc', TRUE, 0),
	c('prop.disc.genes', TRUE, 0),
	c('normal.w.cnas', NA, 0),
	c('normal.cnas', NA, 0),
	c('total.cnas', NA, 0),
	c('prop.disc.genes.os', FALSE, 1),
	c('conc.mean.os', TRUE, 1),
	c('mean.f1score.os', TRUE, 1), 
	c('mean.f1score.gain.os', TRUE, 1), 
	c('mean.f1score.loss.os', TRUE, 1),
	c('mean.ari.smp.os', TRUE, 1),
# c('mean.ari', TRUE, 1),
	c('mean.ari.gene.os', TRUE, 1),
# c('sd.ari', NA, 1)
	c('sd.ari.smp.os', NA, 1),
	c('sd.ari.gene.os', NA, 1)
	);
score.and.sort <- as.data.frame(score.and.sort);
colnames(score.and.sort) <- c('score', 'sort', 'oncoscan');

# state which parameters and scoring metrics are not used for given dataset
skip.params <- c();
if (grepl('nsncnv', proj.stem)) {
	skip.params <- c(skip.params, 'matched');
	};
if (proj.stem == 'bristow' | proj.stem == 'nsncnv_col') {
	skip.params <- c(skip.params, 'col');
	}

skip.scores <- c(
	'mean.ari.gene.os',
# 'sd.ari',
	'mean.f1score.os', #'mean.f1score.gain.os', 'mean.f1score.loss.os',
	'normal.cnas', 'normal.w.cnas', 'total.cnas',
	'sd.inv', 'sd.hk', 'sd.ari.smp.os', 'sd.ari.gene.os',
	'prop.disc.genes'
	);

parameters <- parameters[!parameters %in% skip.params];
score.and.sort <- score.and.sort[!score.and.sort[,1] %in% skip.scores,];

# set colours
none.colour <- 'grey80';
colour.list <- list();
colour.list[['bc']]   	 <- c(none.colour, rev(default.colours(3, palette.type = "spiral.sunrise")));
colour.list[['ccn']]  	 <- c(none.colour, rev(default.colours(2, palette.type = "spiral.afternoon")));
colour.list[['scc']] 	 <- c(none.colour, rev(default.colours(4, palette.type = "div")), "firebrick4");
colour.list[['oth']] 	 <- c(none.colour, rev(default.colours(3, palette.type = "spiral.dusk")));
colour.list[['matched']] <- c('#ABD9E9', '#2C7BB6');
colour.list[['perchip']] <- c(none.colour, 'black');
colour.list[['cnas']] 	 <- c("#f1eef6", "#d4b9da", "#c994c7", "#df65b0", "#dd1c77", "#980043");
colour.list[['col']] 	 <- c(none.colour, 'chartreuse4');

for (i in names(colour.list)) {
	if (!i %in% parameters) {
		colour.list[i] <- NULL;
		}
	}

### FUNCTIONS #################################################################################
# added parameter 'perchip' (not including 'matched' due to too few replicates)

### Load data
# modified from ~/svn/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R
# not sure 'ensemble' applies to current dataset
load.data  <- function(
	dir.name,
	dates,
	total.runs,
	patterns = c('global_*', 'perchip_*'),
	scores = score.and.sort,
	params = parameters,
	use.oncoscan = TRUE
	){

	if (total.runs < 2) { stop("Need at least two samples or results will be wonky!"); }

	if (use.oncoscan) {
		oncoscan.date <- dates[length(dates)];
		dates <- dates[-length(dates)];
		}

	n.params <- length(parameters);
	n.scores <- length(which(score.and.sort[, 3] == 0));
	os.n.scores <- length(which(score.and.sort[, 3] == 1));

	results <- matrix(nrow = total.runs, ncol = (n.scores + os.n.scores));
	params 	<- matrix(nrow = total.runs, ncol = n.params);
	genes	<- list(length = total.runs);
	file.n 	<- 1;

	for (p in 1:length(patterns)) {
		result.patterns <- list.files(
			pattern = paste0('^', patterns[p]),
			include.dirs = TRUE,
			path = dir.name,
			full.names = TRUE
			);

		for (run in 1:length(result.patterns)) {
			if (any(file.exists(paste0(result.patterns[run], '/', dates, '_summary_statistics.txt')))) {
				the.date <- which(file.exists(paste0(result.patterns[run], '/', dates, '_summary_statistics.txt')));
				cur.results <- read.table(
					paste0(result.patterns[run], '/', dates[the.date], '_summary_statistics.txt'),
					header = FALSE
					);
				cur.results <- cur.results[,1:2];
			
				# read in OncoScan comparison results, if available
				if (use.oncoscan & file.exists(paste0(result.patterns[run], '/', oncoscan.date, '_oncoscan_comparison.txt'))) {
					oncoscan.results <- read.table(
						paste0(result.patterns[run], '/', oncoscan.date, '_oncoscan_comparison.txt'),
						header = FALSE
						);
					cur.results <- rbind(cur.results, oncoscan.results, stringsAsFactors = FALSE);
					}

				params[file.n, ]  <- cur.results[match(parameters, cur.results[,2]), 1];
				results[file.n, ] <- cur.results[match(score.and.sort[,1], cur.results[,2]), 1];#  cur.results[, 2] %in% score.and.sort[, 1], 1];
				genes[[file.n]]   <- read.delim(paste0(result.patterns[run], '/', dates[the.date], '_tmr2ref_rounded_counts.txt'));
			} else {
				print(paste("Missing file for", patterns[p], result.patterns[run]));
				}
			file.n <- file.n + 1;				
			}
		}

	last.results <- cur.results;
	
	colnames(params)  <- parameters;
	colnames(results) <- score.and.sort[, 1];

	return(list(
		params = as.data.frame(params),
		scores = as.data.frame(results),
		genes = genes
		));
	}

### Set up covariates for plotting
make.covs <- function(param.combns, border.col = NULL) {
	if (is.null(border.col)) border.col <- 'transparent';

	params.factored <- param.combns;
	colours.factored    <- list();

	for (n in names(param.combns)) {
		levels <- as.numeric(levels(factor(param.combns[, n])));
		params.factored[, n]  <- factor(param.combns[, n], levels = levels);
		colours.factored[[n]] <- colour.list[[n]][levels + 1];
		}

	run.covs <- generate.covariates(
		x = params.factored,
		colour.list = colours.factored,
		col.set = border.col
		);

	return(run.covs);
	}

### Plotting functions
make.pve.barplot <- function(pve.df, fname, var.type = 'Percent Variance Explained') {
	create.barplot(
		ind ~ pve,
		data = pve.df,
		filename = fname,
		sample.order = 'increasing',
		xlab.label = var.type,
		yaxis.lab = pve.df$predictor,
		yaxis.rot = 0,
		xaxis.cex = 1.1,
		yaxis.cex = 1.1,
		xlab.cex = 1.5,
		ylab.cex = 2,
		plot.horizontal = TRUE,
		resolution = 500
		);
	}

make.p.heatmap <- function(pvals, fname) {
	plot.colours <- default.colours(5, palette.type = 'div')[c(1,3,5)];
	plot.labels <- c('p < 0.01', 'p > 0.01', 'p > 0.05');

	pvals.ternary <- pvals;
	pvals.ternary[pvals.ternary >  0.05] <- 3;
	pvals.ternary[pvals.ternary >  0.01 & pvals.ternary <= 0.05] <- 2;
	pvals.ternary[pvals.ternary <= 0.01] <- 1;

	pvals.levels <- levels(as.factor(unlist(pvals.ternary)));
	pvals.levels <- as.numeric(pvals.levels[pvals.levels != "NaN"]);

	pvals.diff <- max(pvals.levels) - min(pvals.levels);
	if (pvals.diff == 1) {
		minimum <- min(pvals.ternary, na.rm = TRUE);
		labels.at <- c(minimum + 0.25, minimum + 0.75);
	} else if (pvals.diff == 2) {
		labels.at <- c(1.33, 2, 2.66);
		}

	create.heatmap(
		x = t(pvals.ternary),
		filename = fname,
		clustering.method = 'none',
		yaxis.lab = rownames(pvals),
		xaxis.lab.top = colnames(pvals), 
		xaxis.cex = 1,
		xaxis.rot = 0,
		yaxis.cex = 1,
		grid.row = TRUE,
		grid.col = TRUE,
		grid.colour = 'gray',
		colour.scheme = plot.colours[min(pvals.levels):max(pvals.levels)],
		total.colours = length(plot.colours[min(pvals.levels):max(pvals.levels)]) + 1,
		colourkey.labels.at = labels.at,
		colourkey.labels = plot.labels[pvals.levels],
		colourkey.cex = 1,
		x.alternating = 3,
		resolution = 500
		);
	}

### Linear model analysis for which params are relevant to run ranking-- evaluate which pipeline parameters are most important
run.glm <- function(glm.data, stem.name, proj.stem) {
	### Generalized linear modelling
	if (grepl('parameters', stem.name)) {
		param.names <- names(glm.data)[-length(names(glm.data))];
		glm.formula <- "log10(rank.prod) ~ ";
		glm.formula <- paste0(glm.formula, paste(param.names, collapse = " + "));
		
		while (length(param.names) > 1) {
			p <- param.names[1];
			param.names <- param.names[-1];
			glm.formula <- paste(glm.formula, "+", paste(p, param.names, sep = "*", collapse = " + "));
			}

		glm.full <- glm(
			formula = glm.formula,
			data = glm.data
			);
	} else {
		metric.names <- names(glm.data)[-length(names(glm.data))];
		glm.formula  <- "log10(rank.prod) ~ ";
		glm.formula  <- paste0(glm.formula, paste(metric.names, collapse = " + "));

		glm.full <- glm(
			formula = glm.formula,
			data = glm.data
			);
		}

	print(glm.formula);

	### Get the percent of variance explained by each parameter
	# pve = percent variance explained from ANOVA performed on glm model
	pve.full 	 <- get.pve(glm.full);
	pve.full$pve <- pve.full$pve * 100;
	pve.full$ind <- seq(1:nrow(pve.full));
	
	make.pve.barplot(pve.full, paste0(plot.dir, generate.filename(stem.name, 'glm_full_pve_barplot', 'png')));
	
	### Model selection by AIC (Akaike information criterion)
	# According to Wikipedia, AIC "is a measure of relative quality of statistical models for a
	# given set of data". It does so by estimating the amount of information lost when a given
	# model is used.
	glm.reduced <- step(glm.full, direction = 'backward');
	pve 		<- get.pve(glm.reduced);
	pve$pve 	<- pve$pve * 100;
	pve$ind 	<- seq(1:nrow(pve));
	
	make.pve.barplot(pve, paste0(plot.dir, generate.filename(stem.name, 'glm_bwelim_pve_barplot', 'png')));

	pdf(file = paste0(plot.dir, generate.filename(stem.name, 'glm_plots', 'pdf')));
	plot(glm.reduced);
	dev.off();

	plot.df <- data.frame(
		resids = stdres(glm.reduced),
		fitted.vals = glm.reduced$fitted.values
		);

	create.scatterplot(
		resids ~ fitted.vals,
		data = plot.df,
		filename = paste0(plot.dir, generate.filename(stem.name, 'glm_bwelim_resid_vs_fitted', 'png')),
		xlab.label = 'fitted values',
		ylab.label = 'Standardized residuals',
		xlimits = c(min(plot.df$fitted.vals) * 0.9, max(plot.df$fitted.vals) * 1.05),
		xlab.cex = 2,
		ylab.cex = 2,
		abline.h = 0,
		abline.col = 'red',
		abline.lty = 2
		);

	# model reduction as shown in http://www.stat.columbia.edu/~martin/W2024/R11.pdf
	glm.reduced2 <- glm(log10(rank.prod) ~ 1, data = glm.data);
	anova(glm.reduced2, glm.full, test = "Chisq");

	return(glm.reduced);
	}

### Random forest analysis in regression mode-- evaluate which scores are most important
run.rf <- function(glm.data, stem.name) {
	if (grepl("parameters", stem.name)) {
		param.names <- names(glm.data)[-length(names(glm.data))];
		glm.formula <- "log10(rank.prod) ~ ";
		glm.formula <- paste0(glm.formula, paste(param.names, collapse = " + "));

		rf.params <- randomForest(
			formula = formula(glm.formula),
			data = glm.data,
			keep.forest = TRUE
			);
	} else {
		metric.names <- names(glm.data)[-length(names(glm.data))];
		glm.formula  <- "log10(rank.prod) ~ ";
		glm.formula  <- paste0(glm.formula, paste(metric.names, collapse = " + "));

		rf.params <- randomForest(
			formula = formula(glm.formula),
			data = glm.data,
			keep.forest = TRUE
			);	# *** add any other metrics here! ***
		}

	param.importance <- importance(rf.params);

	# format for plotting
	param.importance 		   <- as.data.frame(param.importance);
	param.importance$predictor <- rownames(param.importance);
	param.importance$pve 	   <- param.importance$IncNodePurity;
	param.importance$ind 	   <- 1:nrow(param.importance);

	make.pve.barplot(
		pve.df = param.importance,
		fname = paste0(plot.dir, generate.filename(stem.name, 'rf_importance_barplot', 'png')),
		var.type = 'Increase Node Purity'
		);

	return(param.importance);
	}

# switching from default algorithm (Hartigan-Wong) to Lloyd (because inexplicable warning messages..)
cluster.param <- function(param, n.groups, decr) {
	clustered.param <- param;
	clusters <- kmeans(
		na.exclude(param),
		centers = n.groups,
		nstart = 1e3,
		iter.max = 1e4,
		algorithm = "Lloyd"
		)$cluster;

	cluster.order <- order(
		unlist(
			lapply(
				1:n.groups,
				function(k) median(param[clusters==k])
				)
			), decreasing = decr
		);

	for(i in 1:n.groups){
		clustered.param[clusters == cluster.order[i]] <- n.groups^i;
		}

	return(clustered.param);
	}

### MAIN ########################################################################################
# set up directories
main.dir <- ("/.mounts/labs/boutroslab/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV/");

if (grepl('nsncnv', proj.stem)) {
	# if (proj.stem == 'nsncnv_col') dates <- c('2017-01-20', '2017-03-30');
	if (proj.stem == 'nsncnv_col') dates <- c('2017-03-29', '2017-03-31');
	if (proj.stem == 'nsncnv') dates <- c('2017-03-29', '2017-03-31');
	total.runs <- 12672;

	data.dir <- paste0(main.dir, "normalization_assessment/");

	plot.dir <- paste0(main.dir, "plots/total_comparison/");
	if (!file.exists(plot.dir)) dir.create(plot.dir);
	plot.dir <- paste0(plot.dir, paste(dates, collapse = "_"));
	if (!file.exists(plot.dir)) dir.create(plot.dir);
} else if (proj.stem == 'bristow') {
	dates <- c('2017-03-28', '2017-03-31');
	total.runs <- 12672;

	data.dir <- paste0(main.dir, "bristow_assessment/");

	plot.dir <- paste0(main.dir, "bristow_plots/total_comparison/");
	if (!file.exists(plot.dir)) dir.create(plot.dir);
	plot.dir <- paste0(plot.dir, paste(dates, collapse = "_"));
	if (!file.exists(plot.dir))	dir.create(plot.dir);
	}

# collect results data
setwd(data.dir);

# load.data will read in the score files from each directory ** probably very specific to my data **
results <- load.data(
	dir.name = data.dir,
	dates = dates,
	total.runs = total.runs
	);

# order param results by colour list names
results$params <- results$params[, match(names(colour.list), names(results$params))];

# identify and remove any duplicate param combination runs
if (proj.stem == 'bristow') {
	remove.params <- which(apply(results$params, 1, function(x) all(is.na(x))));
} else if (proj.stem == 'nsncnv') {
	remove.params <- which(is.na(results$scores$ari.pts) | is.na(results$scores$prop.disc.genes.os));
} else if (proj.stem == 'nsncnv_col') {
	remove.params <- which(is.na(results$scores$prop.disc.genes.os));
	}

if (length(remove.params) > 0) {
	results$params <- results$params[-remove.params,];
	results$scores <- results$scores[-remove.params,];	
	}

### removing runs where mean.f1score.gain or mean.f1score.loss is NA  ***BE CAREFUL***
if (remove.any.na.runs) {
	to.remove <- which(apply(results$scores, 1, function(x) any(is.na(x))));
	if (length(to.remove) > 0) {
		results$scores <- results$scores[-to.remove,];
		results$params <- results$params[-to.remove,];
		}
	}

{
	# #temp
	# results$scores <- results$scores[which(results$params$cnas < 4),];
	# results$params <- results$params[which(results$params$cnas < 4),];
	# rownames(results$scores) <- rownames(results$params) <- NULL;

	# ### Run k-means on params with many ties *** can decide to re-introduce if needed ***
	# sd.clusters <- cluster.param(
	# 	results$scores$sd.inv.and.hk,
	# 	n.groups = 4,
	# 	decr = FALSE
	# 	);	# have NAs here due to lm so skip this one

	# ari.pts.normcor.clusters <- cluster.param(
	# 	results$scores$ari.pts.normcor,
	# 	n.groups = 8,
	# 	decr = TRUE
	# 	);

	# ari.chip.clusters <- cluster.param(
	# 	results$scores$ari.chip,
	# 	n.groups = 6,
	# 	decr = FALSE
	# 	);

	# if (!'ari.pts.clusters' %in% names(results$scores)) {
	# 	ari.pts.clusters <- cluster.param(
	# 		results$scores$ari.pts,
	# 		n.groups = 8,
	# 		decr = TRUE
	# 		);

	# 	results$scores <- cbind(
	# 		results$scores,
	# 		ari.pts.clusters
	# 		);
	# 	}
}

### Evalute which parameters are associated with each criteria
kw.out  <- matrix(nrow = ncol(results$scores), ncol = ncol(results$params));
aov.out <- matrix(nrow = ncol(results$scores), ncol = ncol(results$params));

for (criteria in 1:ncol(results$scores)) {
	if (colnames(results$scores)[criteria] == 'sd.inv.and.hk') {
		next;
		}

	if (nrow(unique(na.omit(as.vector(results$scores[criteria])))) < 2) {
		print(paste0("Not unique: ", names(results$scores)[criteria]));
		next;
		}

	for (param in 1:ncol(results$params)) {
		if (nrow(unique(na.omit(as.vector(results$params[param])))) > 1) {
			temp.df <- data.frame(
				y = as.vector(t(results$scores[criteria])),
				x = as.vector(t(results$params[param]))
				);
			tryCatch(
				{ kw.out[criteria, param] <- kruskal.test(y ~ x, data = temp.df)$p.value },
				error = function(e) {
					print(paste0(
						"Cannot test at criteria: ", names(results$scores)[criteria], 
						" and param: ", names(results$params)[param]
						));
					kw.out[criteria, param] <- 1;
					}
				); 
			aov.temp <- aov(y ~ x, data = temp.df);
			aov.out[criteria, param] <- summary(aov.temp)[[1]][1, "Pr(>F)"];
		} else {
			kw.out[criteria, param]  <- 1;
			aov.out[criteria, param] <- 1;
			}
		}
	}

kw.out <- as.data.frame(kw.out);
aov.out <- as.data.frame(aov.out);
rownames(kw.out) <- rownames(aov.out) <- colnames(results$scores);
colnames(kw.out) <- colnames(aov.out) <- colnames(results$params);

# Make a heatmap of results
setwd(plot.dir);
if (proj.stem == 'nsncnv_col') {
	make.p.heatmap(aov.out, generate.filename('criteria_x_params_colOnly', 'aov_p_heatmap', 'png'));
	make.p.heatmap(kw.out,  generate.filename('criteria_x_params_colOnly', 'kw_p_heatmap',  'png'));
} else {
	make.p.heatmap(aov.out, generate.filename('criteria_x_params', 'aov_p_heatmap', 'png'));
	make.p.heatmap(kw.out,  generate.filename('criteria_x_params', 'kw_p_heatmap',  'png'));
	}

### Rank data according to each metric
# determine ranks for each variable
setwd(data.dir);

ranks <- results$scores;
sort.decr <- as.vector(score.and.sort$sort[match(names(ranks), as.vector(score.and.sort$score))]);

for (r in 1:ncol(ranks)) {
	if (sort.decr[r]) {
		ranks[, r] <- rank(-1 * results$scores[, r], ties = 'min', na.last = 'keep'); #/ nrow(results$scores);
	} else {
		ranks[, r] <- rank(	  	results$scores[, r], ties = 'min', na.last = 'keep');# / nrow(results$scores);
		}
	}

### Specify 'important' criteria
# DS: originally, EL selected these by taking the top 2 criteria for 'increasing node purity'
if (grepl('nsncnv', proj.stem)) {
	# imp.vars <- c('ari.chip', 'sd.inv.and.hk', 'conc.mean.os', 'mean.ari.smp.os', 'mean.f1score.gain.os');
	imp.vars <- c('ari.pts', 'conc.mean.os', 'mean.ari.smp.os', 'replicates.conc', 'mean.f1score.loss.os', 'mean.f1score.gain.os');
} else if (proj.stem == 'bristow') {
	# imp.vars <- c('mean.f1score.loss.os', 'ari.pts', 'ari.pts.normcor', 'mean.ari.smp.os', 'mean.f1score.gain.os');
	imp.vars <- c('ari.pts.normcor', 'ari.pts', 'replicates.conc', 'conc.mean.os', 'mean.f1score.loss.os');
	}
# imp.vars <- c('replicates.conc', 'ari.pts');
# imp.vars <- c('replicates.conc', 'ari.pts.normcor.clusters');
# imp.vars <- c('replicates.conc', 'prop.disc.genes', 'ari.pts');
# imp.vars <- qw('ari.chip ari.type');
# imp.vars <- qw('conc.mean.os prop.disc.genes.os');
# imp.vars <- qw('ari.chip ari.pts.normcor ari.type replicates.conc ari.pts ari.pts.clusters conc.mean.os prop.disc.genes.os');
# imp.vars <- qw('ari.chip ari.pts.normcor ari.type replicates.conc ari.pts');

# DS: when matched param is not excluded (see above) all top ranking runs correspond to matched = 1 and ari.pts is NA 
# geometric mean
# used when running glm and rf below
overall.rank <- apply(
	X = ranks,
	MARGIN = 1,
	FUN = function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); }
	);

overall.rank <- apply(
	X = ranks[, imp.vars, drop = FALSE],
	MARGIN = 1,
	FUN = function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); }
	);

ranks$overall <- rank(overall.rank);
run.orders    <- order(ranks$overall);

### Repeat analysis with each variable omitted and track ranks of runs
reduced.ranks <- ranks[, -tail(1:ncol(ranks), n = 1)];

for (r in 1:ncol(reduced.ranks)) {
	reduced.ranks[,r] <- rank(
		apply(
			X = ranks[, -r], 
			MARGIN = 1,
			FUN = function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); }
			),
		ties = 'min',
		na.last = 'keep'
		);
	}

# summarize # top ranks
final.rank 		<- apply(reduced.ranks, 1, function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); } );
n.top.ranks 	<- apply(ranks, 		1, function(f) { length(which(f == 1)); } );
n.top5.ranks 	<- apply(ranks, 		1, function(f) { length(which(f <= 5)); } );
# n.top5imp.ranks <- apply(ranks[, imp.vars, drop = FALSE], 1, function(f) { length(which(f <= 5)); } );# this doesn't get accessed anywhere

# add the rank prod to the matrix
reduced.ranks$rank.prod <- final.rank;

### Track various score for runs
ordering.scores.df <- data.frame(
	top  = n.top.ranks * -1,	# need -1 for rand prod below
	top5 = n.top5.ranks * -1,
	imp.rank.prod = overall.rank,
	# imp.rank.prod = apply(ranks[, imp.vars, drop = FALSE], 1, function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); } ),# same as overall.rank
	all.rank.prod = apply(ranks[, !colnames(ranks) %in% 'overall'], 1, function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); } ),
	# all.rank.prod = apply(ranks, 1, function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); } ),
	reduced.ranks.score = final.rank
	);

reduced.ranks$all.rank.prod <- ordering.scores.df$all.rank.prod;
reduced.ranks$imp.rank.prod <- ordering.scores.df$imp.rank.prod;

ordering.scores.df$summary.rank.prod <- apply(
	X = apply(ordering.scores.df, 2, rank),
	MARGIN = 1,
	FUN = function(f) { prod(f) ^ (1 / length(f)); }
	);

### reorder runs based on meta-ranks
run.orders2 <- order(
	ordering.scores.df$imp.rank.prod,
	ordering.scores.df$all.rank.prod,
	ordering.scores.df$top,
	ordering.scores.df$top5
	);

reduced.ranks <- reduced.ranks[run.orders2,];
reduced.ranks$overall <- seq(1:nrow(reduced.ranks));

ordering.scores.df <- ordering.scores.df[run.orders2,];
ordering.scores.df$meta.rank <- seq(1:nrow(ordering.scores.df));

# put variables back to positive space
ordering.scores.df$top  <- ordering.scores.df$top  * -1;
ordering.scores.df$top5 <- ordering.scores.df$top5 * -1;

### Output some scores to file
if (proj.stem == 'nsncnv_col') {
	fname <- generate.filename('params_x_scores_colOnly', 'ranked_output', 'txt')
} else {
	fname <- generate.filename('params_x_scores', 'ranked_output', 'txt')
	}

write.table(
	cbind(results$params[run.orders2,], ordering.scores.df, ranks[run.orders2,]),
	file = fname,
	sep= "\t",
	quote = FALSE,
	row.names = FALSE
	);

### Fit some linear models to evaluate statistically which variables are important for the rank product
params.stem <- 'parameters';
criteria.stem <- 'criteria';

if (proj.stem == 'nsncnv_col') {
	params.stem <- 'parameters_colOnly';
	criteria.stem <- 'criteria_colOnly';
	}

# for parameters
glm.df <- data.frame(lapply(results$params, factor));
glm.df$rank.prod <- final.rank;
glm.params <- run.glm(glm.df, stem.name = params.stem, proj.stem = proj.stem);
# rf.params  <- run.rf(glm.df,  stem.name = 'parameters', proj.stem = proj.stem);

# glm.df.unmatched <- data.frame(lapply(results$params[-which(is.na(results$scores$normal.cnas)),], factor));
# glm.df.unmatched$rank.prod <- final.rank[-which(is.na(results$scores$normal.cnas))];
# glm.params.unmatched <- run.glm(glm.df.unmatched, stem.name = 'parameters_unmatched');
# # rf.params.unmatched  <- run.rf(glm.df.unmatched,  stem.name = 'parameters_unmatched');

# for ranking criteria
glm.df.criteria <- results$scores;
glm.df.criteria$rank.prod <- overall.rank;
# glm.df.criteria$rank.prod <- overall.rank.all;# not just 'imp.vars'
# glm.criteria <- run.glm(glm.df.criteria, stem.name = 'criteria');

if (!remove.any.na.runs & proj.stem == 'nsncnv_col') {
	glm.df.criteria <- glm.df.criteria[apply(glm.df.criteria, 1, function(x) !any(is.na(x))),];
	}

if (any(apply(glm.df.criteria, 1, function(x) any(is.na(x))))) {
	stop("Problem with glm.df.criteria!");
	}

rf.criteria <- run.rf(glm.df.criteria, stem.name = criteria.stem);# won't run with missing values

# glm.df.criteria.unmatched <- results$scores[-which(is.na(results$scores$normal.cnas)),];
# glm.df.criteria.unmatched$rank.prod <- overall.rank[-which(is.na(results$scores$normal.cnas))];
# # glm.criteria.unmatched <- run.glm(glm.df.criteria.unmatched, stem.name = 'criteria_unmatched');
# rf.criteria.unmatched  <- run.rf(glm.df.criteria.unmatched, stem.name = 'criteria_unmatched');

# {
# 	glm.df.binomial <- apply(glm.df.unmatched,2,as.numeric);
# 	cols.to.transform <- which(!colnames(glm.df.binomial) %in% c('cnas', 'rank.prod'));
# 	glm.df.binomial[, cols.to.transform][glm.df.binomial[, cols.to.transform] > 0] <- 1;
# 	glm.df.binomial <- as.data.frame(glm.df.binomial);
# 	glm.df.binomial[,-ncol(glm.df.binomial)] <- apply(glm.df.binomial[,-ncol(glm.df.binomial)],2,as.factor);

# 	glm.binomial <- glm(
# 		log10(rank.prod)/10 ~ perchip + bc + ccn + scc + oth + cnas + col + 
# 			perchip*bc + perchip*ccn + perchip*scc + perchip*oth + perchip*cnas + perchip*col +
# 			bc*ccn + bc*scc + bc*oth + bc*cnas + bc*col + 
# 			ccn*scc + ccn*oth + ccn*cnas + ccn*col + 
# 			scc*oth + scc*cnas + scc*col + 
# 			oth*cnas + oth*col + 
# 			cnas*col,
# 		data = glm.df.binomial,
# 		family = 'binomial'
# 		);

# 	pve.binomial 	 <- get.pve(glm.binomial);
# 	pve.binomial$pve <- pve.full$pve * 100;
# 	pve.binomial$ind <- seq(1:nrow(pve.binomial));
		
# 	glm.binomial.reduced 	 <- step(glm.binomial, direction = 'backward');
# 	pve.binomial.reduced 	 <- get.pve(glm.reduced);
# 	pve.binomial.reduced$pve <- pve$pve * 100;
# 	pve.binomial.reduced$ind <- seq(1:nrow(pve));

# 	pdf(file = paste0(plot.dir, generate.filename(stem.name, 'glm_plots', 'pdf')));
# 	plot(glm.binomial.reduced);
# 	dev.off();

# 	plot.df <- data.frame(
# 		resids = stdres(glm.binomial.reduced),
# 		fitted.vals = glm.binomial.reduced$fitted.values
# 		);

# 	create.scatterplot(
# 		resids ~ fitted.vals,
# 		data = plot.df,
# 		# filename = paste0(plot.dir, generate.filename(stem.name, 'glm_bwelim_resid_vs_fitted', 'png')),
# 		xlab.label = 'fitted values',
# 		ylab.label = 'Standardized residuals',
# 		xlimits = c(min(plot.df$fitted.vals) * 0.9, max(plot.df$fitted.vals) * 1.05),
# 		xlab.cex = 2,
# 		ylab.cex = 2,
# 		abline.h = 0,
# 		abline.col = 'red',
# 		abline.lty = 2
# 		);

# 	# model reduction as shown in http://www.stat.columbia.edu/~martin/W2024/R11.pdf
# 	glm.binomial.reduced2 <- glm(log10(rank.prod) ~ 1, data = glm.df.binomial);
# 	anova(glm.binomial.reduced2, glm.binomial, test = "Chisq");
# }

### Plotting ################################################################
setwd(plot.dir);

# make a plot of run scores showing params as covariates
# create legend
border.col <- 'black';
run.legend <- list(
	col = list(
		colours = rev(colour.list[['col']]),
		title = 'Collapsed Probeset',
		border = border.col,
		cex = 1,
		labels = c('yes', 'no')
		),
	cnas = list(
		colours = rev(colour.list[['cnas']]),
		title = 'CNA Method',
		border = border.col,
		cex = 1,
		labels = c('KD: 0.996, 0.695, 0.73, 0.956', 'KD: 0.983, 0.748, 0.727, 0.969', 'KD: 0.998, 0.79, 0.88, 0.989', 'KD (default): 0.85, 0.95', 'Normal Min/Max Thresholds', 'NS-inspired Thresholds')
		),
	oth = list(
		colours = rev(colour.list[['oth']]),
		title = 'Other Norm',
		border = border.col,
		cex = 1,
		labels = c('quantile', 'rank.normal', 'vsn', 'none')
		),
	matched = list(
		colours = rev(colour.list[['matched']]),
		title = 'Reference',
		border = border.col,
		cex = 1,
		labels = c('matched', 'pooled')
		),
	scc = list(
		colours = rev(colour.list[['scc']]),
		title = 'Sample Content',
		border = border.col,
		cex = 1,
		labels = c('invariant.norm', 'low.cv.geo.mean', 'top.geo.mean', 'total.sum', 'housekeeping.geo.mean', 'none')
		),
	bc = list(
		colours = rev(colour.list[['bc']]),
		title = 'Background',
		border = border.col,
		labels.cex = 1,
		labels = c('max', 'mean.2sd', 'mean', 'none')
		),
	ccn = list(
		colours = rev(colour.list[['ccn']]),
		title = 'Code Count',
		border = border.col,
		labels.cex = 0.75,
		labels = c('geo.mean', 'sum', 'none')
		),
	perchip = list(
		colours = rev(colour.list[['perchip']]),
		title = 'Per Chip',
		border = border.col,
		labels.cex = 1,
		labels = c('yes', 'no')
		)
	);

# remove parameters from legend that don't apply to proj.stem
for (n in names(run.legend)) {
	if (!n %in% names(colour.list)) {
		run.legend[[n]] <- NULL;
		} 
	}
names(run.legend) <- rep("legend", length(run.legend));

# heatmap of ranks
if (proj.stem == 'nsncnv_col') {
	fname <- generate.filename('runs_summary_colsOnly', 'ranks_heatmap', 'png');
} else {
	fname <- generate.filename('runs_summary', 'ranks_heatmap', 'png');
	}

ranks.plotting <- ranks[run.orders, ];
run.covs <- make.covs(results$params[run.orders,]);
create.heatmap(
	x = ranks.plotting,
	filename = fname,
	yaxis.lab = TRUE,
	yaxis.cex = 1,
	xaxis.lab = NULL,
	cluster.dimensions = 'none',
	fill.colour = none.colour,
	covariates.top = run.covs,
	covariate.legends = run.legend,
	legend.title.cex = 0.75,
	legend.cex = 0.6,
	legend.side = 'right',
	colour.scheme = c('purple', 'white'),
	colourkey.labels.at = c(1, 3000, 6000, 9000, 12000),
	colourkey.cex = 1,
	resolution = 600,
	height = 7
	);

# order matrix according to scores
scores.plotting <- as.data.frame(results$scores);
scores.plotting$run.orders <- run.orders;

for (j in names(scores.plotting)) {
	if (j == 'run.orders') {
		ordering <- order(scores.plotting[,j]);
	} else {
		ordering <- order(scores.plotting[,j], decreasing = score.and.sort$sort[score.and.sort$score %in% j]);
		}
	scores.plotting <- scores.plotting[ordering,];
	run.covs <- make.covs(results$params[ordering,]);

	scores.plotting[scores.plotting < 0] <- 0;

	for (i in 1:ncol(scores.plotting)) {
		scores.plotting[,i] <- scores.plotting[,i]/max(scores.plotting[,i], na.rm = TRUE);
		}

	if (proj.stem == 'nsncnv_col') {
		fname <- generate.filename(paste0('runs_summary_', j, '_colsOnly'), 'heatmap', 'png');
	} else {
		fname <- generate.filename(paste0('runs_summary_', j), 'heatmap', 'png');
		}

	create.heatmap(
		x = scores.plotting,
		filename = fname,
		yaxis.lab = TRUE,
		yaxis.cex = 1,
		xaxis.lab = NULL,
		clustering.method = 'none',
		cluster.dimensions = 'none',
		fill.colour = none.colour,
		covariates.top = run.covs,
		scale.data = TRUE,
		covariate.legends = run.legend,
		legend.title.cex = 0.85,
		legend.cex = 0.6,
		legend.side = 'right',
		colour.scheme = c('yellow', 'white', 'purple'),
		colourkey.cex = 0.6,
		resolution = 600,
		height = 9,
		width = 10
		);
	}

### For each criteria, plot the distribution in the cohort
final.params <- results$params[run.orders2,];
final.scores <- results$scores[run.orders2,];

run.covs2 <- make.covs(final.params);

key.labels <- rbind(
	c('col', 'Collapsed'),
	c('cnas', 'CNAmeth'),
	c('oth', 'OthNorm'),
	c('scc', 'SmpCont'),
	c('bc', 'BkgdCor'),
	c('ccn', 'CodeCnt'),
	c('perchip', 'PerChip'),
	c('matched', 'Reference')
	);

for (x in 1:ncol(final.scores)) {
	temp.params <- final.params[order(final.scores[, x], decreasing = TRUE),];
	temp.scores <- final.scores[order(final.scores[, x], decreasing = TRUE),];
	temp.covs2  <- make.covs(temp.params);

	temp.df <- data.frame(
		ranks = 1:nrow(temp.scores),
		value = temp.scores[, x]
		);

	if (proj.stem == 'nsncnv_col') {
		fname <- generate.filename(colnames(final.scores)[x], 'colsOnly_distribution', 'png');
	} else {
		fname <- generate.filename(colnames(final.scores)[x], 'distribution', 'png')
		}

	create.barplot(
		formula = value ~ ranks,
		data = temp.df,
		filename = fname,
		ylab.label = colnames(final.scores)[x],
		ylab.cex = 2,
		xaxis.tck = 0,
		xaxis.lab = rep('', nrow(final.scores)),
		col = 'black',
		border.col = 'black',
		legend = list(
			bottom = list(
				fun = covariates.grob(
					covariates = rev(temp.covs2),
					side = 'top',
					size = 0.6,
					ord = seq(1:nrow(temp.df))
					)
				),
			right = list(
				fun = legend.grob(
					legends = run.legend,
					title.just = 'left',
					title.cex = .75,
					label.cex = .75,
					size = 2
					)
				)
			),
		key = list(
			x = 1,
			y = -0.037,
			text = list(
				lab = rev(paste(
					key.labels[match(names(temp.params), key.labels[,1]), 2],
					signif(kw.out[x, names(temp.params)], digits = 2),
					sep = ": "
					)),
				cex = 0.66
				),
			padding.text = 1
			),
		bottom.padding = 7,
		resolution = 600,
		height = 8.5,
		);
	}

### Make a multiplot for reduced rank analysis
### heatmap of reduced ranks
hmap <- create.heatmap(
	x = reduced.ranks,
#	generate.filename('runs_summary', 'reduced_ranks_heatmap', 'png'),
	yaxis.lab = TRUE,
	xaxis.lab = rep('', nrow(reduced.ranks)),
	cluster.dimensions = 'none',
	# covariates.top = run.covs2,
	# covariate.legends = run.legend,
	# legend.side = 'right',
	# legend.cex = 0.6,
	colour.scheme = c('purple', 'white'),
	# colourkey.cex = 1,
	# resolution = 1200,
	# height = 9
	);

### Get summary data to use in barplots
top.ranks.bplot <- create.barplot(
	top ~ meta.rank,
	data = ordering.scores.df,
#	filename = generate.filename('reduced_ranks', 'n_top_ranks_barplot', 'png'),
	ylab.label = "# top ranks",
	main = "# top ranks",
	lwd = 0,
	#xlimits = c(0.5,nrow(reduced.ranks)+ 0.5),
	xlimits = c(0, nrow(reduced.ranks)),
	xaxis.tck = 0,
	ylimits = c(0, max(ordering.scores.df$top, na.rm = TRUE) + 0.5),
	yat = pretty(c(0, max(ordering.scores.df$top, na.rm = TRUE))),
	# yat = seq(0, 20, by = 3),
	xat = rep('', nrow(ordering.scores.df))
	);

top5.ranks.bplot <- create.barplot(
	top5 ~ meta.rank,
	data = ordering.scores.df,
#	filename = generate.filename('reduced_ranks', 'n_top5_ranks_barplot', 'png'),
	lwd = 0,
	ylab.label = "# top 5 ranks",
	main = "# top 5 ranks",
	ylimits = c(0, max(ordering.scores.df$top, na.rm = TRUE) + 0.5),
	yat = pretty(c(0, max(ordering.scores.df$top5, na.rm = TRUE))),
	# yat = seq(0, 20, by = 3),
	xlimits = c(0,nrow(reduced.ranks)),
	#xlimits = c(0.5,nrow(reduced.ranks) + 0.5),
	xaxis.tck = 0,
	xat = rep('', nrow(ordering.scores.df))
	);

# try to evaluate criteria
rank.dist <- vector(length = ncol(reduced.ranks));
for(criteria in 1:ncol(reduced.ranks)){
	# distance between ranks without criteria and final ranks
	rank.dist[criteria] <- dist(t(reduced.ranks[, c(criteria, ncol(reduced.ranks))]));
	}

rank.dist.df <- data.frame(
	criteria = seq(1:length(rank.dist)),
	distance = rank.dist
	);

criteria.barplot <- create.barplot(
	criteria ~ distance,
	data = rank.dist.df,
#	filename = generate.filename("reduced_ranks", "criteria_distance_to_rank_barplot", 'png'),
	plot.horizontal = TRUE,
	lwd = 0,
	yaxis.lab = rep('', nrow(rank.dist.df)),
	xlab.label = 'Euclidean distance',
	main = 'Euclidean distance',
	xaxis.lab = c('0', "1\n\t\tEuclidean Distance (100 000s)", "2", '3'),
	xat = seq(0, 300000, 100000),
	xlimits = c(0, 310000),
	ylimits = c(0.5, nrow(rank.dist.df) + 0.5)
	);

# evaluate whether there is a bias in parameters between top quartile and remaining runs
n.runs <- nrow(final.params);
top.quartile <- floor(n.runs / 3);
prop.out <- matrix(nrow = ncol(final.params), ncol = 3);

for (p in 1:ncol(final.params)) {
	# split data based on number of levels for current parameter
	top.x <- floor(n.runs / length(unique(final.params[, p])));
	prop.out[p,] <- as.numeric(prop.test(
		x = cbind(
			table(factor(
				final.params[1:(top.x), p],
				levels = t(unique(final.params[p])))),
			table(factor(
				final.params[seq(from = (top.x + 1), to = n.runs), p],
				levels = t(unique(final.params[p]))
				))
			)
		)[qw("statistic parameter p.value")]);
	}

colnames(prop.out) <- qw("x2 df p");
prop.out <- as.data.frame(prop.out);
prop.out$log10p <- -log10(prop.out$p);
prop.out$run <- seq(1:nrow(prop.out));

prop.out$log10p[prop.out$log10p == Inf] <- max(prop.out$log10p[prop.out$log10p != Inf], na.rm = TRUE) * 1.1;

params.barplot <- create.barplot(
	run ~ log10p,
	data = prop.out,
#	filename = generate.filename("reduced_ranks", "parameter_prop_p_barplot", 'png'),
	plot.horizontal = TRUE,
	as.table = TRUE,
	yaxis.lab = rep('', nrow(prop.out)),
	yat = seq(1, nrow(prop.out)),
	lwd = 0,
	xlab.label = '-log10(p)',
	main = '-log10(p)',
	abline.v = -log10(0.05),
	abline.col = 'red',
	xlimits = c(0, 21),
	# xlimits = c(0, max(prop.out$log10p)),
	# xat = seq(0, max(prop.out$log10p), by = 5),
	xaxis.lab = c('0', '5', '10\n\t   -log10(p)', '15', '20+'),
	xaxis.tck = 0.5,
	yaxis.tck = 0,
	ylimits = c(0.5, nrow(prop.out) + 0.5)
	);

### make a heatmap of run parameters
cov.colours <- unique(as.vector(unlist(colour.list)));

run.params.matrix <- final.params;
if ('matched' %in% parameters) run.params.matrix$matched <- run.params.matrix$matched + 1;
if ('cnas' %in% parameters) run.params.matrix$cnas <- run.params.matrix$cnas + 1;

max.val <- max(run.params.matrix[, 1], na.rm = TRUE);
for (param in 2:ncol(run.params.matrix)) {
	cur.vals <- run.params.matrix[, param];
	cur.vals[which(cur.vals > 0)] <- cur.vals[which(cur.vals > 0)] + max.val;
	max.val <- max(c(max.val, cur.vals), na.rm = TRUE);
	run.params.matrix[, param] <- cur.vals;
	}

covs.hmap <- create.heatmap(
	x = run.params.matrix,
#	filename = generate.filename('reduced_ranks', 'covariates', 'png'),
	total.colours = max(run.params.matrix, na.rm = TRUE) + 2,
	colour.scheme = cov.colours[1:(max(run.params.matrix, na.rm = TRUE) + 1)],
	clustering.method = 'none',
	xaxis.tck = 0,
	xaxis.lab = rep('', nrow(run.params.matrix)),
	yaxis.tck = 0,
	yaxis.cex = 0.8,
	grid.col = FALSE,
	grid.row = TRUE,
	row.colour = 'gray', 
	col.colour = 'gray',
	force.grid.row = TRUE,
	print.colour.key = FALSE,
	resolution = 1200
	);

### Put all together
if (proj.stem == 'nsncnv_col') {
	fname <- generate.filename("reduced_ranks_colsOnly", "multiplot", "png");
} else {
	fname <- generate.filename("reduced_ranks", "multiplot", "png");
	}

create.multiplot(
	plot.objects = list(top5.ranks.bplot, top.ranks.bplot, hmap, criteria.barplot, covs.hmap, params.barplot),
	filename = fname,
	panel.heights = c(2.3,6,1,1),
	panel.widths = c(3,1),
	plot.layout = c(2, 4),
	xaxis.tck = c(0),
	# xaxis.tck = c(0.5, 0,0,0,0,0),
	yaxis.tck = 0.5,
	y.spacing = c(-1,-1,0),
	x.spacing = -0.5,
	layout.skip = c(F,T,F,T,F,F,F,F),
	retrieve.plot.labels = TRUE,
	left.padding = 7,
	xlimits = list(
		c(1,nrow(reduced.ranks)),
		c(1,nrow(reduced.ranks)),
		c(1,nrow(reduced.ranks)),
		c(0, 100),
		c(1,nrow(reduced.ranks)),
		c(0,nrow(ranks))
		),
	y.relation = 'free',
	xaxis.cex = 0.55,
	# xlab.label = c('', '-log10(p)', '', '    Euclidean distance', '', ''),
	ylab.label = c("\t\t\t", "\t\t\t", "\t\t\t", 'Top 5 ranks   Top ranks       '),
	ylab.padding = 1,
	ylab.cex = 0.65,
	bottom.padding = 1,
	xlab.cex = 1,
	yaxis.cex = 0.6,
	x.relation = 'sliced',
	legend = run.legend,
	resolution = 500,
	width = 6.5
	);

### Make a heatmap showing the ranks of the top methods
top.n <- 50;
all.scores <- cbind(ordering.scores.df[1:top.n,], ranks[run.orders2[1:top.n] , ]);
all.scores[, c(1:2)] <- apply(all.scores[, c(1:2)], 2, function(f) 1 - (f / max(f)));# need to transform top and top5 rankings
all.scores[,-c(1:2)] <- apply(all.scores[,-c(1:2)], 2, function(f) f / max(f));
if(top.n < 200) {
	top.run.covs <- make.covs(results$params[run.orders2[1:top.n],], border.col = 'black');
} else {
	top.run.covs <- make.covs(results$params[run.orders2[1:top.n],], border.col = NULL);	
	}

# all.scores <- apply(all.scores, 2, function(x) { if (all(x == 0)) { x <- rep(1, nrow(all.scores)); } else { x; } });
# all.scores <- apply(all.scores, 2, function(x) { if (all(x == 1)) { x <- rep(min(all.scores), nrow(all.scores)); } else { x; } });

### Put all together
if (proj.stem == 'nsncnv_col') {
	fname <- generate.filename('top_runs_colsOnly', 'metric_comparison', 'png');
} else {
	fname <- generate.filename('top_runs', 'metric_comparison', 'png');
	}

create.heatmap(
	x = t(all.scores),
	filename = fname,
	cluster.dimensions = 'none',
	xaxis.lab = colnames(all.scores),
	yaxis.lab = 1:top.n,
	xaxis.cex = 1,
	yaxis.cex = 0.7,
	covariates = top.run.covs,
	covariate.legends = run.legend,
	scale.data = TRUE,
	legend.between.row = 0.5,
	legend.border.padding = 0.5,
	legend.cex = 0.5,
	legend.title.cex = 0.6,
	colour.scheme = c('white', 'black'),
	colour.alpha = 0.5,
	colourkey.labels = c('Best', paste0(nrow(all.scores), 'th')),
	colourkey.labels.at = c(0, max(all.scores)),
	colourkey.cex = 1.5,
	total.colours = 99,
	resolution = 600,
	top.padding = 15,
	height = 10
	);

### SESSION_INFO ##################################################################################
if (proj.stem == 'nsncnv_col') {
	save.session.profile(paste0(plot.dir, generate.filename('Run_Comparison_colOnly', 'Session_Info','txt')), FALSE);
} else {
	save.session.profile(paste0(plot.dir, generate.filename('Run_Comparison', 'Session_Info','txt')), FALSE);
	}
