### compare_preprocessing_results_2.0.R ###########################################################

### PREAMBLE ######################################################################################
library(mclust);
library(NanoStringNorm);
library(BoutrosLab.utilities);
library(BoutrosLab.pipeline.limma);
library(BoutrosLab.plotting.general);
#library(BoutrosLab.statistics.general);
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

# set up directories
main.dir <- ("/.mounts/labs/boutroslab/private/AlgorithmEvaluations/microarrays/NanoStringNormCNV/");
data.dir <- paste0(main.dir, "normalization_assessment/");
plot.dir <- paste0(main.dir, "plots/");

# set colours
colour.list <- list();
none.colour <- 'grey80';
colour.list[['bc']]   	 <- c(none.colour, rev(default.colours(3, palette.type = "spiral.sunrise")));
colour.list[['ccn']]  	 <- c(none.colour, rev(default.colours(2, palette.type = "spiral.afternoon")));
colour.list[['scc']] 	 <- c(none.colour, rev(default.colours(5, palette.type = "div")));
colour.list[['oth']] 	 <- c(none.colour, rev(default.colours(3, palette.type = "spiral.dusk")));
colour.list[['matched']] <- c('#ABD9E9', '#2C7BB6');
colour.list[['perchip']] <- c(none.colour, 'black');
colour.list[['cnas']] 	 <- c("#f1eef6", "#d4b9da", "#c994c7", "#df65b0", "#dd1c77", "#980043");
# colour.list[['cnas']] 	 <- c(default.colours(5, palette.type = "spiral.dawn")[5:2]);
colour.list[['col']] 	 <- c(none.colour, 'chartreuse4');

### FUNCTIONS #################################################################################
# added parameters 'perchip'

### Load data
# modified from ~/svn/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R
# not sure 'ensemble' applies to current dataset
load.data  <- function(
	dir.name,
	dates,
	total.runs,
	patterns = c('global_*', 'perchip_*'),
	n.scores = 15,
	n.params = 9
	){

	if (total.runs < 2) { stop("Need at least two samples or results will be wonky!"); }

	results <- matrix(nrow = total.runs, ncol = (n.scores - n.params));
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
			if (file.exists(paste0(result.patterns[run], '/', dates, '_summary_statistics.txt'))) {
				cur.results <- read.table(
					paste0(result.patterns[run], '/', dates, '_summary_statistics.txt'),
					header = FALSE
					);

				results[file.n, ] <- as.numeric(t(cur.results[(n.params + 1):n.scores, 1]));
				params[file.n, ]  <- t(cur.results[1:n.params, 1]);
				genes[[file.n]]   <- read.delim(paste0(result.patterns[run], '/', dates, '_tmr2ref_rounded_counts.txt'));
			} else {
				print(paste("Missing file for", patterns[p], result.patterns[run]));
				}
			file.n <- file.n + 1;				
			}
		}

	last.results <- cur.results;
	
	colnames(params)  <- last.results[1:n.params, 2];
	colnames(results) <- last.results[(n.params + 1):n.scores, 2];

	cols.to.remove <- qw("sd.inv sd.hk cand.gene.cor prop.disc.genes lmyc.validation normals.w.cnas");
	# cols.to.remove <- qw("sd.inv sd.hk ari.type ari.pts.normcor cand.gene.cor prop.disc.genes lmyc.validation");

	results <- results[, -which(colnames(results) %in% cols.to.remove)];

	return(list(
		params = as.data.frame(params),
		scores = as.data.frame(results),
		genes = genes
		));
	}

### Set up covariates for plotting
make.covs <- function(parameters) {
	run.covs <- generate.covariates(
		x = data.frame(
		 	perchip = factor(parameters[, "perchip"], levels = c(0, 1)),
			ccn 	= factor(parameters[, "ccn"], 	  levels = c(0, 1, 2)),
			bc 		= factor(parameters[, "bc"], 	  levels = c(0, 1, 2, 3)),
			scc 	= factor(parameters[, "scc"], 	  levels = c(0, 1, 2, 3, 4, 5)),
			matched = factor(parameters[, 'matched'], levels = c(0, 1)),
			oth 	= factor(parameters[, 'oth'], 	  levels = c(0, 1, 2, 3)),
			cnas 	= factor(parameters[, 'cnas'], 	  levels = as.numeric(levels(as.factor(parameters$cnas)))),
			col 	= factor(parameters[, 'col'], 	  levels = c(0, 1))
			),
		colour.list = list(
	 		perchip = colour.list[['perchip']],
			ccn 	= colour.list[['ccn']],
			bc 		= colour.list[['bc']],
			scc 	= colour.list[['scc']],
			matched = colour.list[['matched']],
			oth 	= colour.list[['oth']],
			cnas 	= colour.list[['cnas']][1:nlevels(as.factor(results$params$cnas))],
			col 	= colour.list[['col']]
			)#,
		# col.set = 'black'
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
run.glm <- function(glm.data, stem.name) {
	### Generalized linear modelling
	if (stem.name == 'parameters') {
		glm.full <- glm(
			log10(rank.prod) ~ perchip + bc + ccn + scc + matched + oth + cnas + col + 
				perchip*bc + perchip*ccn + perchip*scc + perchip*matched + perchip*oth + perchip*cnas + perchip*col +
				bc*ccn + bc*scc + bc*matched + bc*oth + bc*cnas + bc*col + 
				ccn*scc + ccn*matched + ccn*oth + ccn*cnas + ccn*col + 
				scc*matched + scc*oth + scc*cnas + scc*col + 
				matched*oth + matched*cnas + matched*col + 
				oth*cnas + oth*col + 
				cnas*col,
			data = glm.data
			);
	} else {
		glm.full <- glm(
			log10(rank.prod) ~ ari.chip + ari.pts.normcor + ari.type + ari.pts + sd.inv.and.hk + replicates.conc,
			data = glm.data
			);
		}

	### Get the percent of variance explained by each parameter
	pve.full 	 <- get.pve(glm.full);
	pve.full$pve <- pve.full$pve * 100;
	pve.full$ind <- seq(1:nrow(pve.full));
	
	make.pve.barplot(pve.full, generate.filename(stem.name, 'glm_full_pve_barplot', 'png'));
	
	### Model selection by AIC (Akaike information criterion)
	# According to Wikipedia, AIC "is a measure of relative quality of statistical models for a
	# given set of data". It does so by estimating the amount of information lost when a given
	# model is used.
	glm.reduced <- step(glm.full, direction = 'backward');
	pve 		<- get.pve(glm.reduced);
	pve$pve 	<- pve$pve * 100;
	pve$ind 	<- seq(1:nrow(pve));
	
	make.pve.barplot(
		pve,
		paste0(plot.dir, generate.filename(stem.name, 'glm_bwelim_pve_barplot', 'png'))
		);
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
		xlab.cex = 2,
		ylab.cex = 2
		);

	return(glm.reduced);
	}

### Random forest analysis in regression mode-- evaluate which scores are most important
run.rf <- function(glm.data, stem.name) {
	if (stem.name == 'parameters') {
		rf.params <- randomForest(
			formula = log10(rank.prod) ~ bc + ccn + scc + oth + matched + perchip + cnas + col,
			# formula = log10(rank.prod) ~ bc + ccn + scc + matched + oth + cnas + col,
			data = glm.data,
			keep.forest = TRUE
			);
	} else {
		rf.params <- randomForest(
			formula = log10(rank.prod) ~ ari.chip + ari.pts.normcor + ari.type + ari.pts + sd.inv.and.hk + replicates.conc,
			# formula = log10(rank.prod) ~ ari.chip + ari.pts + sd.inv.and.hk + replicates.conc,
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
		param.importance,
		paste0(plot.dir, generate.filename(stem.name, 'rf_importance_barplot', 'png')),
		var.type = 'Increase Node Purity'
		);

	return(param.importance);
	}

### MAIN ########################################################################################
# collect results data
setwd(data.dir);

# load.data will read in the score files from each directory ** probably very specific to my data **
results <- load.data(
	dir.name = data.dir,
	dates = '2017-01-18',
	total.runs = 12672,
	n.scores = 20,
	n.params = 8
	);

# order param results by colour list names
results$params <- results$params[, match(names(colour.list), names(results$params))];

{
	# colnames(results$params)[8] <- 'collapsed';
	# print(results$scores[, grep(x = colnames(results$scores), pattern = 'validation')]);

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

	# results$scores <- cbind(
	# 	results$scores,
	# 	cand.gene.clusters,
	# 	replicates.conc.clusters,
	# 	total.val.clusters
	# 	);
}

### Evalute which parameters are associated with each criteria
kw.out  <- matrix(nrow = ncol(results$scores), ncol = ncol(results$params));
aov.out <- matrix(nrow = ncol(results$scores), ncol = ncol(results$params));

for (criteria in 1:ncol(results$scores)) {
	if (colnames(results$scores)[criteria] == 'sd.inv.and.hk') {
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
					kw.out[criteria, param] <- NA;
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
rownames(kw.out) <- colnames(results$scores);
colnames(kw.out) <- colnames(results$params);

aov.out <- as.data.frame(aov.out);
rownames(aov.out) <- colnames(results$scores);
colnames(aov.out) <- colnames(results$params);

# Make a heatmap of results
setwd(plot.dir);
make.p.heatmap(aov.out, generate.filename('criteria_x_params', 'aov_p_heatmap', 'png'));
make.p.heatmap(kw.out,  generate.filename('criteria_x_params', 'kw_p_heatmap',  'png'));

### Rank data according to each metric
# determine ranks for each variable
setwd(data.dir);
# ranks <- results$scores[, -ncol(results$scores)];
ranks <- results$scores[, !(colnames(results$scores) %in% c('normal.cnas', 'total.cnas'))];# this modification may be incorrect

# sort.decr <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE);#EL
sort.decr <- c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE);#DS --see svn/Collaborators/RobBristow/nanostring_validation/normalization/compare_preprocessing_results.R r35400
print(cbind(metric = names(ranks), decreasing = sort.decr));

for (r in 1:ncol(ranks)) {
	if (sort.decr[r]) {
		ranks[, r] <- rank(-1 * results$scores[, r], ties = 'min', na.last = 'keep'); #/ nrow(results$scores);
	} else {
		ranks[, r] <- rank(	  	results$scores[, r], ties = 'min', na.last = 'keep');# / nrow(results$scores);
		}
	}

### Specify 'important' criteria
# imp.vars <- c('replicates.conc', 'ari.pts.normcor.clusters');
imp.vars <- c('replicates.conc', 'ari.pts');
# imp.vars <- c('ari.chip', 'ari.pts.normcor', 'ari.type', 'sd.inv.and.hk', 'replicates.conc', 'ari.pts');

overall.rank <- apply(
	X = ranks[, imp.vars, drop = FALSE],
	MARGIN = 1,
	FUN = function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); }
	);

ranks$overall <- rank(overall.rank);
run.orders    <- order(ranks$overall);

### Repeat analysis with each variable omitted and track ranks of runs
reduced.ranks <- ranks[, -ncol(ranks)];
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
n.top5imp.ranks <- apply(ranks[, imp.vars, drop = FALSE], 1, function(f) { length(which(f <= 5)); } );# this doesn't get accessed anywhere

# add the rank prod to the matrix
reduced.ranks$rank.prod <- final.rank;

### Track various score for runs
ordering.scores.df <- data.frame(
	top  = n.top.ranks * -1,	# need -1 for rand prod below
	top5 = n.top5.ranks * -1,
	imp.rank.prod = apply(ranks[, imp.vars, drop = FALSE], 1, function(f) { prod(f) ^ (1 / length(f)); } ),# same as overall.rank
	all.rank.prod = apply(ranks, 1, function(f) { prod(na.omit(f)) ^ (1 / length(na.omit(f))); } ),
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

reduced.ranks 		  <- reduced.ranks[run.orders2,];
reduced.ranks$overall <- seq(1:nrow(reduced.ranks));

ordering.scores.df 	   <- ordering.scores.df[run.orders2,];
ordering.scores.df$run <- seq(1:nrow(ordering.scores.df));

reduced.ranks$my.rank <- rank(apply(reduced.ranks[, -ncol(reduced.ranks)], 1, sum), ties = 'average');
ordering.scores.df$my.rank <- rank(apply(ordering.scores.df[, -ncol(ordering.scores.df)], 1, sum), ties = 'average');

# put variables back to positive space
ordering.scores.df$top  <- ordering.scores.df$top  * -1;
ordering.scores.df$top5 <- ordering.scores.df$top5 * -1;

### Output some scores to file
write.table(
	cbind(results$params[run.orders2,], ordering.scores.df, ranks[run.orders2,]),
	generate.filename('params_x_scores', 'ranked_output', 'txt'),
	sep= "\t",
	quote = FALSE,
	row.names = FALSE
	);

### Fit some linear models to evaluate statistically which variables are important for the rank product
# for parameters
glm.df 			 <- data.frame(lapply(results$params, factor));
glm.df$rank.prod <- final.rank;
glm.params 		 <- run.glm(glm.df, stem.name = 'parameters');
rf.params 		 <- run.rf(glm.df,  stem.name = 'parameters');

# for ranking criteria
glm.df.criteria 		  <- results$scores;
glm.df.criteria$rank.prod <- overall.rank;
glm.criteria 	  		  <- run.glm(glm.df.criteria, stem.name = 'criteria');
rf.criteria 	 		  <- run.rf(glm.df.criteria,  stem.name = 'criteria');

### Plotting ################################################################
setwd(plot.dir);

# make a plot of run scores showing params as covariates
# order matrix according to scores
scores.plotting <- as.data.frame(results$scores[run.orders, ]);
ranks.plotting  <- ranks[run.orders, ];

run.covs <- make.covs(results$params[run.orders,]);

border.col <- 'black';
run.legend <- list(
	legend = list(
		colours = rev(colour.list[['col']]),
		title = 'Collapsed Probeset',
		border = border.col,
		cex = 1,
		labels = c('yes', 'no')
		),
	legend = list(
		colours = rev(colour.list[['cnas']]),
		title = 'CNA Method',
		border = border.col,
		cex = 1,
		labels = c('KD: 0.89, 0.69, 0.65, 0.87', 'KD: 0.98, 0.84, 0.92, 0.97', 'KD: 0.998, 0.79, 0.88, 0.989', 'KD (default): 0.85, 0.95', 'Normal Min/Max Thresholds', 'NS-inspired Thresholds')
		),
	legend = list(
		colours = rev(colour.list[['oth']]),
		title = 'Other Norm',
		border = border.col,
		cex = 1,
		labels = c('quantile', 'rank.normal', 'vsn', 'none')
		),
	legend = list(
		colours = rev(colour.list[['matched']]),
		title = 'Reference',
		border = border.col,
		cex = 1,
		labels = c('matched', 'pooled')
		),
	legend = list(
		colours = rev(colour.list[['scc']]),
		title = 'Sample Content',
		border = border.col,
		cex = 1,
		labels = c('invariant.norm', 'low.cv.geo.mean', 'top.geo.mean', 'total.sum', 'housekeeping.geo.mean', 'none')
		),
	legend = list(
		colours = rev(colour.list[['bc']]),
		title = 'Background',
		border = border.col,
		labels.cex = 1,
		labels = c('max', 'mean.2sd', 'mean', 'none')
		),
	legend = list(
		colours = rev(colour.list[['ccn']]),
		title = 'Code Count',
		border = border.col,
		labels.cex = 0.75,
		labels = c('geo.mean', 'sum', 'none')
		),
	legend = list(
		colours = rev(colour.list[['perchip']]),
		title = 'Per Chip',
		border = border.col,
		labels.cex = 1,
		labels = c('yes', 'no')
		)
	);

create.heatmap(
	x = scores.plotting[,! names(scores.plotting) %in% c('normal.cnas', 'total.cnas')],
	generate.filename('runs_summary', 'heatmap', 'png'),
	yaxis.lab = TRUE,
	yaxis.cex = 1,
	xaxis.lab = NULL,
	cluster.dimensions = 'none',
	fill.colour = none.colour,
	covariates.top = run.covs,
	scale.data = TRUE,
	covariate.legends = run.legend,
	legend.title.cex = 0.85,
	legend.cex = 0.75,
	legend.side = 'right',
	colour.scheme = c('purple', 'white', 'yellow'),
	resolution = 1000,
	height = 9,
	width = 10
	);

### heatmap of ranks
create.heatmap(
	x = ranks.plotting,
	generate.filename('runs_summary', 'ranks_heatmap', 'png'),
	yaxis.lab = TRUE,
	yaxis.cex = 1,
	xaxis.lab = NULL,
	cluster.dimensions = 'none',
	fill.colour = none.colour,
	covariates.top = run.covs,
	covariate.legends = run.legend,
	legend.title.cex = 0.75,
	legend.cex = 0.75,
	legend.side = 'right',
	colour.scheme = c('white', 'purple'),
	colourkey.cex = 1,
	resolution = 600,
	height = 9
	);

### For each criteria, plot the distribution in the cohort
final.params <- results$params[run.orders2,];
final.scores <- results$scores[run.orders2,];

run.covs2 <- make.covs(final.params);

for (x in 1:ncol(final.scores)) {
	temp.params <- final.params[order(final.scores[, x], decreasing = TRUE),];
	temp.scores <- final.scores[order(final.scores[, x], decreasing = TRUE),];
	temp.covs2  <- make.covs(temp.params);

	temp.df <- data.frame(
		ranks = 1:nrow(temp.scores),
		value = temp.scores[, x]
		);

	create.barplot(
		formula = value ~ ranks,
		data = temp.df,
		filename = generate.filename(colnames(final.scores)[x], 'distribution', 'png'),
		ylab.label = colnames(final.scores)[x],
		ylab.cex = 2,
		xaxis.tck = 0,
		xaxis.lab = rep('', nrow(final.scores)),
		col = 'black',# bar.colours,
		border.col = 'black',#bar.colours,
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
					# size = 1.5,
					title.cex = .75,
					label.cex = .75,
					size = 2
					)
				)
			),
		key = list(
			x = 1,
			y = -0.017,
			text = list(
				lab = paste(
					c('Collapsed', 'CNAmeth', 'OthNorm', 'Reference', 'SmpCont', 'BkgdCor', 'CodeCnt', 'PerChip'),
					signif(kw.out[x, c('col', 'cnas', 'oth', 'matched', 'scc', 'bc', 'ccn', 'perchip')], digits = 2),
					sep = ": "
					),
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
	colour.scheme = c('white', 'purple'),
	# colourkey.cex = 1,
	# resolution = 1200,
	# height = 9
	);

### Get summary data to use in barplots
top.ranks.bplot <- create.barplot(
	top ~ run,
	data = ordering.scores.df,
#	filename = generate.filename('reduced_ranks', 'n_top_ranks_barplot', 'png'),
	ylab.label = "# top ranks",
	main = "# top ranks",
	lwd = 0,
	#xlimits = c(0.5,nrow(reduced.ranks)+ 0.5),
	xlimits = c(0, nrow(reduced.ranks)),
	xaxis.tck = 0,
#	ylimits = c(0, 6),
	yat = pretty(c(0, max(ordering.scores.df$top, na.rm = TRUE))),
	# yat = seq(0, 20, by = 3),
	xat = rep('', nrow(ordering.scores.df))
	);

top5.ranks.bplot <- create.barplot(
	top5 ~ run,
	data = ordering.scores.df,
#	filename = generate.filename('reduced_ranks', 'n_top5_ranks_barplot', 'png'),
	lwd = 0,
	ylab.label = "# top 5 ranks",
	main = "# top 5 ranks",
#	ylimits = c(0, 12),
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
	xaxis.lab = c('0', "3", "6\nEuclidean Distance (100 000s)", '9', '12'),
	xat = seq(0, 1200000, 300000),
	xlimits = c(0, max(rank.dist.df$distance)),
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
	xlimits = c(0, max(prop.out$log10p)),
	xat = seq(0, max(prop.out$log10p), by = 5),
	xaxis.lab = c('0', '5', '10', '15\n\t   -log10(p)', '20', '25', '30', '35'),
	xaxis.tck = 0.5,
	yaxis.tck = 0,
	ylimits = c(0.5, nrow(prop.out) + 0.5)
	);

### make a heatmap of run parameters
cov.colours <- unique(as.vector(unlist(colour.list)));

run.params.matrix <- final.params;
run.params.matrix$matched <- run.params.matrix$matched + 1;
run.params.matrix$cnas <- run.params.matrix$cnas + 1;

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
create.multiplot(
	plot.objects = list(top5.ranks.bplot, top.ranks.bplot, hmap, criteria.barplot, covs.hmap, params.barplot),
	filename = generate.filename("reduced_ranks", "multiplot", "png"),
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
	# ylimits = list(
	# 	c(1.5, 4.5),
	# 	NULL,
	# 	c(0.5, ncol(reduced.ranks) + 0.5),
	# 	NULL,
	# 	c(0.5, 4.5),
	# 	c(0.5, 4.5)
	# 	),
ylimits = list(NULL,NULL,NULL,NULL,NULL,c(.5, 9.5)),
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
all.scores <- cbind(ordering.scores.df[1:10 , ], ranks[run.orders2[1:10] , ]);
all.scores[, c(1:2)] <- apply(all.scores[, c(1:2)], 2, function(f) 1 - (f / max(f)));	# need to transform top and top5 rankings
all.scores[, -c(1:2)] <- apply(all.scores[, -c(1:2)], 2, function(f) f / max(f));
top.run.covs <- make.covs(results$params[run.orders2[1:10],]);

### Put all together
create.heatmap(
	x = t(all.scores),
	filename = generate.filename('top_runs', 'metric_comparison', 'png'),
	cluster.dimensions = 'none',
	xaxis.lab = colnames(all.scores),
	yaxis.lab = 1:10,
	xaxis.cex = 1,
	yaxis.cex = 1,
	covariates = top.run.covs,
	covariate.legends = run.legend,
	colour.scheme = c('white', 'black'),
	colour.alpha = 0.5,
	colourkey.labels = c('Best', '10th'),
	colourkey.labels.at = c(0.01, 1),
	colourkey.cex = 1.5,
	total.colours = 99,
	resolution = 600
	);
