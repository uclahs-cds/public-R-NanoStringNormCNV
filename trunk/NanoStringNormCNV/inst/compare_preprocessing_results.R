
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
load_all("~/blab_svn/BoutrosLab/Resources/code/R/prostate.acgh.biomarkers");
load_all("~/blab_svn/BoutrosLab/Resources/code/R/NanoStringNormCNV/");
source("~/blab_svn/BoutrosLab/Training/elalonde/OncoScan_reprocess/cna.plotting.functions.R");
source("~/blab_svn/BoutrosLab/Resources/code/R/ParameterEval/R/generate.covariates.R")
source("~/blab_svn/BoutrosLab/Resources/code/R/BoutrosLab.statistics.general/R/get.pve.R")
source("~/svn/BoutrosLab/Collaborators/RobBristow/nanostring_validation/normalization/accessory_functions.R")


### FUNCTIONS #################################################################################
make.covs <- function(parameters){
	run.covs <- generate.covariates(
		x = data.frame(
			bc = factor(parameters[, "bc"], levels = c(0, 1)),
			cc = factor(parameters[, "cc"], levels = c(0, 1)),
			scc = factor(parameters[, "scc"], levels = c(1,2,3,4)),
			ref = factor(parameters[, 'matched'], levels = c(0, 1)),
			other = factor(parameters[, 'other'], levels = c(0,1,2,3)),
			cnas = factor(parameters[, 'cnas'], levels = c(0,1,2,3,4))
			collapsed = factor(parameters[, 'collapsed'], levels = c(0,1))
			),
		colour.list = list(
			bc = c('white','black'),
			cc = c('white','blue'),
			scc = c('yellow', 'orange', 'red', 'darkmagenta'),
			ref = c('white', 'pink'),
			other = brewer.pal(4, 'Set1'),
			cnas = c('white', 'violet','purple', 'darkorchid4', 'slateblue4')
			collapsed = c('white', 'indianred3')
			)
#		col.set = 'white'
		);
	return(run.covs);
	}

### Random forest analysis in regression mode -- evaluate which scores are most important
run.rf <- function(glm.data, stem.name){
	# remove ensemle data since there are many missing values
	glm.data <- glm.data[-nrow(glm.data) , ];
	if(stem.name == 'parameters'){
		rf.params <- randomForest(formula=log10(rank.prod)~bc+cc+scc+matched+other+cnas+collapsed, data = glm.data, keep.forest=TRUE);
	}else{
		rf.params <- randomForest(formula=log10(rank.prod)~ari.chip+ari.pts+sd.inv.and.hk+replicates.conc , data = glm.data, keep.forest=TRUE);	# *** add any other metrics here! ***
		}
	param.importance <- importance(rf.params);
	# format for plotting
	param.importance <- as.data.frame(param.importance);
	param.importance$predictor <- rownames(param.importance);
	param.importance$pve <- param.importance$IncNodePurity;
	param.importance$ind <- 1:nrow(param.importance);
	make.pve.barplot(param.importance, generate.filename(stem.name, 'rf_importance_barplot', 'png'), var.type = 'Increase Node Purity');
	return(param.importance);
	}

### Linear model analysis for which params are relevant to run ranking-- evaluate which pipeline parameters are most important
run.glm <- function(glm.data, glm.formula){
# exclude the other param
	if(length(unique(glm.data$other)) == 1){	# *** these if statements check if multiple options were specified for certain parameters- if not, don't include in linear model
		if(length(unique(glm.data$matched)) == 1){
			if(length(unique(na.omit(glm.data$collapsed))) == 1){
				glm.full <- glm(log10(rank.prod) ~ cc + scc + cnas + bc +
					cc*scc + cc*cnas + cc*bc +
					scc*cnas + scc*bc +
					cnas * bc,
					data = glm.data
					);
			}else{
				glm.full <- glm(log10(rank.prod) ~ bc + cc + scc + cnas + collapsed +
					bc*cc + bc*scc + bc*cnas + bc*collapsed + 
					cc*scc + cc*cnas + cc*collapsed +
					scc*cnas + scc*collapsed +
					cnas * collapsed,
					data = glm.data
					);
				}
		}else{
			glm.full <- glm(log10(rank.prod) ~ bc + cc + scc + matched + cnas + collapsed +
				bc*cc + bc*scc + bc*matched + bc*cnas + bc*collapsed +
				cc*scc + cc*matched + cc*cnas + cc*collapsed +
				scc*matched + scc*cnas + scc*collapsed +
				matched*cnas + matched*collapsed +
				cnas*collapsed,
				data = glm.data
				);
			}
	}else{
		if(length(unique(glm.data$matched)) == 1){
			glm.full <- glm(log10(rank.prod) ~ bc + cc + scc + other + cnas + collapsed +
				bc*cc + bc*scc + bc*other + bc*cnas + bc*collapsed +
				cc*scc + cc*other + cc*cnas + cc*collapsed +
				scc*other + scc*cnas + scc*collapsed +
				other*cnas + other*collapsed +
				cnas*collapsed,
				data = glm.data
				);
		}else{
			glm.full <- glm(log10(rank.prod) ~ bc + cc + scc + matched + other + cnas + collapsed +
				bc*cc + bc*scc + bc*matched + bc*other + bc*cnas + bc*collapsed +
				cc*scc + cc*matched + cc*other + cc*cnas + cc*collapsed +
				scc*matched + scc*other + scc*cnas + scc*collapsed +
				matched*other + matched*cnas + matched*collapsed +
				other*cnas + other*collapsed +
				cnas*collapsed,
				data = glm.data
				);
			}
		}
	pve.full <- get.pve(glm.full);
	pve.full$pve <- pve.full$pve * 100;
	pve.full$ind <- seq(1:nrow(pve.full));
	make.pve.barplot(pve.full, generate.filename(stem.name, 'glm_full_pve_barplot', 'png'));
	
	glm.reduced <- step(glm.full, direction = 'backward');
	pve <- get.pve(glm.reduced);
	pve$pve <- pve$pve * 100;
	pve$ind <- seq(1:nrow(pve));
	make.pve.barplot(pve, generate.filename(stem.name, 'glm_bwelim_pve_barplot', 'png'));
	pdf(file=generate.filename(stem.name, 'glm_plots', 'pdf'));
	plot(glm.reduced);
	dev.off();
	plot.df <- data.frame(resids = stdres(glm.reduced), fitted.vals = glm.reduced$fitted.values);
	create.scatterplot(resids ~ fitted.vals, data = plot.df, filename=generate.filename(stem.name, 'glm_bwelim_resid_vs_fitted', 'png'), xlab.label = 'fitted values', ylab.label = 'Standardized residuals', xlab.cex = 2, ylab.cex = 2);
	return(glm.reduced);
	}

make.pve.barplot <- function(pve.df, fname, var.type = 'Percent Variance Explained'){
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

make.p.heatmap <- function(pvals, fname){
	pvals.ternary <- pvals;
	pvals.ternary[pvals.ternary > 0.05] <- 3;
	pvals.ternary[pvals.ternary <= 0.01] <- 1;
	pvals.ternary[pvals.ternary <= 0.05] <- 2;
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
		colour.scheme = c('red','pink', 'white'),
		total.colours = 4,
		colourkey.labels.at = c(1.33, 2, 2.66),
		colourkey.labels = c('p < 0.01', 'p < 0.05', 'p > 0.05'),
		colourkey.cex = 1,
		x.alternating = 3,
		resolution = 500
		);
	}

### MAIN ########################################################################################
### collect result data
data.path <- ("~/isilon/private/Collaborators/RobBristow/cna_biomarkers/validation/4_Nanostring/data/normalization_assessment_outliers_removed_wlm_wother/analysis");
setwd(data.path);
# load.data will read in the score files from each directory-- probably very specific to my data
results <- load.data(dir.name = paste0(data.path, "/../"), dates = '2015-07-01', total.runs = 490, n.scores = 32, load.ensemble = F);
colnames(results$params)[8] <- 'collapsed';
print(results$scores[, grep(x=colnames(results$scores), pattern = 'validation')]);


### Run k-meams on params with many ties, *** can decide to re-introduce if needed ***
#sd.clusters <- cluster.param(results$scores$sd.inv.and.hk, n.groups = 4, decr = FALSE);	# have NAs here due to lm so skip this one
#ari.pts.normcor.clusters <- cluster.param(results$scores$ari.pts.normcor, n.groups = 8, decr = TRUE);
#ari.chip.clusters <- cluster.param(results$scores$ari.chip, n.groups = 6, decr = FALSE);
#results$scores <- cbind(results$scores, cand.gene.clusters, replicates.conc.clusters, total.val.clusters);


### Evalute which parameters are associated with each criteria
kw.out <- matrix(nrow = ncol(results$scores), ncol = ncol(results$params));
aov.out <- matrix(nrow = ncol(results$scores), ncol = ncol(results$params));
for (criteria in 1:ncol(results$scores)){
	if(colnames(results$scores)[criteria] == 'sd.inv.and.hk'){
		next;
		}
	for (param in 1:ncol(results$params)){
		if(nrow(unique(na.omit(as.vector(results$params[param])))) > 1){
			temp.df <- data.frame(y=as.vector(t(results$scores[criteria])), x = as.vector(t(results$params[param])));
			kw.out[criteria, param] <- kruskal.test(y ~ x, data = temp.df)$p.value;
			aov.temp <- aov(y ~ x, data = temp.df);
			aov.out[criteria, param] <- summary(aov.temp)[[1]][1, "Pr(>F)"];
		}else{
			kw.out[criteria, param] <- 1;
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
make.p.heatmap(aov.out, generate.filename('criteria_x_params', 'aov_p_heatmap', 'png'));
make.p.heatmap(kw.out, generate.filename('criteria_x_params', 'kw_p_heatmap', 'png'));


### Rank data according to each metric
# determine ranks for each variable
ranks <- results$scores[ , -ncol(results$scores)];
print(colnames(ranks));
sort.decr <- c(T,  T, F, T, T, F);
print(sort.decr)
for(r in 1:ncol(ranks)){
	if(sort.decr[r]){
		ranks[, r] <- rank(-1*results$scores[,r], ties = 'min', na.last = 'keep'); #/ nrow(results$scores);
	}else{
		ranks[, r] <- rank(results$scores[,r], ties = 'min', na.last = 'keep');# / nrow(results$scores);
		}
	}

### Specify 'important' criteria
imp.vars <- c('replicates.conc', 'ari.pts.normcor.clusters'); 
overall.rank <- apply(ranks[,imp.vars], 1, function(f) prod(na.omit(f))^(1/length(na.omit(f))));
ranks$overall <- rank(overall.rank);
run.orders <- order(ranks$overall);

### Repeat analysis with each variable omitted and track ranks of runs
reduced.ranks <- ranks[ , -ncol(ranks)];
for(r in 1:ncol(reduced.ranks)){
	reduced.ranks[,r] <- rank(apply(ranks[ , -r], 1, function(f) prod(na.omit(f))^(1/length(na.omit(f)))), ties = 'min', na.last = 'keep');
	}
# summarize # top ranks
final.rank <- apply(reduced.ranks, 1, function(f) prod(na.omit(f))^(1/length(na.omit(f))));
n.top.ranks <- apply(ranks,1,function(f) length(which(f == 1)));
n.top5.ranks <- apply(ranks,1,function(f) length(which(f <= 5)));
n.top5imp.ranks <- apply(ranks[, imp.vars],1,function(f) length(which(f <= 5)));
# add the rank prod to the matrix
reduced.ranks$rank.prod <- final.rank;


### Track various score for runs
ordering.scores.df <- data.frame(
	top = n.top.ranks * -1,	# need -1 for rand prod below
	top5 = n.top5.ranks * -1,
	imp.rank.prod = apply(ranks[, imp.vars], 1, function(f) prod(f)^(1/length(f))),
	all.rank.prod = apply(ranks, 1, function(f) prod(na.omit(f))^(1/length(na.omit(f)))),
	reduced.ranks.score = final.rank
	);
reduced.ranks$all.rank.prod <- ordering.scores.df$all.rank.prod;
reduced.ranks$imp.rank.prod <- ordering.scores.df$imp.rank.prod;
ordering.scores.df$summary.rank.prod <- apply(apply(ordering.scores.df, 2, rank), 1,function(f) prod(f)^(1/length(f)));
### reorder runs based on meta-ranks
run.orders2 <- order(ordering.scores.df$imp.rank.prod, ordering.scores.df$all.rank.prod, ordering.scores.df$top, ordering.scores.df$top5);
reduced.ranks <- reduced.ranks[run.orders2 , ];
reduced.ranks$overall <- seq(1:nrow(reduced.ranks));
ordering.scores.df <- ordering.scores.df[run.orders2 , ];
ordering.scores.df$run <- seq(1:nrow(ordering.scores.df));
# put variables back to positive space
ordering.scores.df$top <- ordering.scores.df$top * -1;
ordering.scores.df$top5 <- ordering.scores.df$top5 * -1;

### Output some scores to file
write.table(cbind(results$params[run.orders2,], ordering.scores.df, ranks[run.orders2,]), generate.filename('params_x_scores', 'ranked_output', 'txt'), sep= "\t", quote = F);

### Fit some linear models to evaluate statistically which variables are important for the rank product
# for parameters
glm.df <- data.frame(lapply(results$params, factor));
glm.df$rank.prod <- final.rank;
rf.params <- run.rf(glm.df, stem.name ='parameters');

# for ranking criteria
glm.df.criteria <- results$scores;
glm.df.criteria$rank.prod <- overall.rank;
glm.reduced.criteria <- run.glm(glm.df.criteria, stem.name = 'criteria');

### Plotting ################################################################
# make a plot of run scores showing params as covariates
# order matrix according to scores
scores.plotting <- as.data.frame(results$scores[run.orders, ]);
ranks.plotting <- ranks[run.orders, ];
run.covs <- make.covs(results$params[run.orders,]);
run.legend <- list(
	legend = list(
		colours = c('black', 'white'),
		title = 'Background',
		border = 'white',
		labels.cex = 1,
		labels = c("mean.2sd", "none")
		),
	legend = list(
		colours = c('blue', 'white'),
		title = 'Code Count',
		border = 'white',
		labels.cex = 0.75,
		labels = c('geo.mean', 'none')
		),
	legend = list(
		colours = c('yellow', 'orange', 'red', 'darkmagenta'),
		title = 'Sample Content',
		border = 'white',
		cex = 1,
		labels = c('housekeeping', 'top.geo.mean', 'invariant', 'inv-lm')
		),
#	legend = list(
#		colours = c('pink', 'white'),
#		title = 'Reference',
#		border = 'white',
#		cex = 1,
#		labels = c('Matched', 'Pooled')
#		),
	legend = list(
		colours = brewer.pal(4, 'Set1'), #c('green', 'white'),
		title = 'Other',
		border = 'white',
		cex = 1,
		labels = c('None', 'VSN', 'Quantile', 'Rank normal')
		),
	legend = list(
		colours = c('white', 'violet', 'purple', 'darkorchid4', 'slateblue4'),
		title = 'CNA method',
		border = 'white',
		cex = 1,
		labels = c('Thresholds','KD-normals', 'KD 0.9,0,87,0.93,0.96', 'KD 0.92,0.87,0.97,0.99', 'KD 0.9,0.885,0.91,0.94')
#		),
#	legend = list(
#		colours = c('white', 'indianred3'),
#		title = 'Collapsed',
#		border = 'white',
#		cex = 1,
#		labels = c('No','Yes')
		)
	);

create.heatmap(
	x = scores.plotting,
	generate.filename('runs_summary', 'heatmap', 'png'),
	yaxis.lab = TRUE,
	xaxis.lab = NULL,
	cluster.dimensions = 'none',
	covariates.top = run.covs,
	scale.data=T,
	covariate.legends = run.legend,
	colour.scheme = c('purple', 'white', 'yellow'),
	resolution = 600
	);

### heatmap of ranks
create.heatmap(
	x = ranks.plotting,
	generate.filename('runs_summary', 'ranks_heatmap', 'png'),
	yaxis.lab = TRUE,
	xaxis.lab = NULL,
	cluster.dimensions = 'none',
	covariates.top = run.covs,
	covariate.legends = run.legend,
	colour.scheme = c('white', 'purple'),
	resolution = 600
	);

### For each criteria, plot the distribution in the cohort
final.params <- results$params[run.orders2 , ];
final.scores <- results$scores[run.orders2 , ];

run.covs2 <- make.covs(final.params);
for(x in 1:ncol(final.scores)){
	temp.params <- final.params[order(final.scores[,x], decreasing = T) , ];
	temp.scores <- final.scores[order(final.scores[,x], decreasing = T) , ];
	temp.covs2 <- make.covs(temp.params);
	temp.df <- data.frame(
		ranks = 1:nrow(temp.scores),
		value = temp.scores[,x]
		);
	#temp.df <- data.frame(
	#	ranks = final.scores$overall,
	#	value = final.scores[,x]
	#	);

	create.barplot(
		formula = value ~ ranks,
#		sample.order = 'decreasing',
		data = temp.df,
		filename = generate.filename(colnames(final.scores)[x], 'distribution', 'png'),
		ylab.label = colnames(final.scores)[x],
		ylab.cex = 2,
		xaxis.tck = 0,
		xaxis.lab = rep('', nrow(final.scores)),
		legend = list(
			bottom = list(fun = covariates.grob(covariates = temp.covs2, side = 'top', size = 0.8, ord = seq(1:nrow(temp.df)))),
#			bottom = list(fun = covariates.grob(covariates = run.covs2, side = 'top', size = 0.8, ord = seq(1:nrow(temp.df)))),
			right = list(fun = legend.grob(legends = run.legend, title.just = 'left'))
			),
		key = list(
			x = 1,
			y =0.002,
			text = list(lab = paste(c('Background', 'Code Count', 'Sample Content', 'Other', 'CNAs'), signif(kw.out[x , ], digits = 2))),
			#text = list(lab = paste(c('Background', 'Code Count', 'Sample Content', 'Reference', 'Other', 'CNAs', 'Collapsed'), signif(kw.out[x , ], digits = 2))),
			padding.text = 1
			),
		bottom.padding = 6,
		resolution = 500
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
	covariates.top = run.covs2,
	covariate.legends = run.legend,
	colour.scheme = c('white', 'purple'),
	resolution = 600
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
	xlimits = c(0,nrow(reduced.ranks)),
	xaxis.tck = 0,
#	ylimits = c(0,6),
	yat = seq(0,20,by=3),
	xat = rep('', nrow(ordering.scores.df))
	);
top5.ranks.bplot <- create.barplot(
	top5 ~ run,
	data = ordering.scores.df,
#	filename = generate.filename('reduced_ranks', 'n_top5_ranks_barplot', 'png'),
	lwd = 0,
	ylab.label = "# top 5 ranks",
	main = "# top 5 ranks",
#	ylimits = c(0,12),
	yat = seq(0,20,by=3),
	xlimits = c(0,nrow(reduced.ranks)),
	#xlimits = c(0.5,nrow(reduced.ranks) + 0.5),
	xaxis.tck = 0,
	xat = rep('', nrow(ordering.scores.df))
	);


# try to evaluate criteria
rank.dist <- vector(length = ncol(reduced.ranks));
for(criteria in 1:ncol(reduced.ranks)){
	# distance between ranks without criteria and final ranks
	rank.dist[criteria] <- dist(t(reduced.ranks[ , c(criteria, ncol(reduced.ranks))]));
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
	xaxis.lab = c('0', '20', "40\nEuclidean Distance  ", "60", '80', '100', '120', '140'),
	xlimits = c(0, max(rank.dist.df$distance) + 5),
	ylimits = c(0.5,nrow(rank.dist.df)+0.5)
	);


# evaluate whether there is a bias in parameters between top quartile and remaining runs
n.runs <- nrow(final.params);
top.quartile <- floor(n.runs / 3);
prop.out <- matrix(nrow = ncol(final.params), ncol = 3);
for(p in 1:ncol(final.params)){
	# split data based on number of levels for current parameter
	top.x <- floor(n.runs / length(unique(final.params[,p])));
	prop.out[p, ] <- as.numeric(prop.test(x = cbind(table(factor(final.params[1:(top.x), p], levels = t(unique(final.params[p])))), table(factor(final.params[seq(from= (top.x+1), to = n.runs), p], levels = t(unique(final.params[p]))))))[qw("statistic parameter p.value")]);
	}
colnames(prop.out) <- qw("x2 df p");
prop.out <- as.data.frame(prop.out);
prop.out$log10p <- -log10(prop.out$p);
prop.out$run <- seq(1:nrow(prop.out));

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
	xlimits = c(0, max(prop.out$log10p) + 0.5),
	xat = seq(0,14,by=2),
	xaxis.lab = c('0', '2\n\t   -log10(p)', '4', '6', '8', '10', '12', '14'),
	xaxis.tck = 0.5,
	ylimits = c(0.5,6.5)
	);

### make a heatmap of run parameters
cov.colours <- qw("white black blue yellow orange red darkmagenta #E41A1C #377EB8 #4DAF4A #984EA3 violet purple darkorchid4 slateblue4");# indianred3");
#cov.colours <- qw("white black blue yellow orange red darkmagenta pink green violet purple darkorchid4");

run.params.matrix <- final.params;
run.params.matrix$other <- run.params.matrix$other+1;
#run.params.matrix$perchip <- ifelse(run.params.matrix$perchip == 1, 0, 1); # reverse this to make colours work
max.val <- 1;
for(param in 2:ncol(run.params.matrix)){
	cur.vals <- run.params.matrix[ , param];
	print(cur.vals);
	cur.vals[which(cur.vals > 0)] <- cur.vals[which(cur.vals > 0)] + max.val;
	max.val <- max(c(max.val, cur.vals), na.rm=T);
	run.params.matrix[, param] <- cur.vals;
	}

covs.hmap <- create.heatmap(
	x = run.params.matrix,
#	filename = generate.filename('reduced_ranks', 'covariates', 'png'),
	total.colours = max(run.params.matrix, na.rm=T) + 2,
	colour.scheme = cov.colours[1:(max(run.params.matrix, na.rm = T)+1)],
	clustering.method = 'none',
	xaxis.tck = 0,
	xaxis.lab = rep('', nrow(run.params.matrix)),
	yaxis.tck = 0,
	yaxis.cex = 0.8,
	grid.col = FALSE,
	grid.row = TRUE,
	row.colour = 'gray', 
	col.colour = 'gray',
	print.colour.key = FALSE,
	resolution = 500
	);


### Put all together
create.multiplot(
	plot.objects = list(top5.ranks.bplot, top.ranks.bplot, hmap, criteria.barplot, covs.hmap, params.barplot),
	filename = generate.filename("reduced_ranks", "multiplot", "png"),
	panel.heights = c(1,6,1,1),
	panel.widths = c(3,1),
	plot.layout = c(2, 4),
	xaxis.tck = 0.5,
	yaxis.tck = 0.5,
	y.spacing = c(-1,-1,0),
	x.spacing = -0.5,
	layout.skip = c(F,T,F,T,F,F,F,F),
	retrieve.plot.labels = TRUE,
	left.padding = 7,
	xlimits = list(c(1,nrow(reduced.ranks)), c(1,nrow(reduced.ranks)), c(1,nrow(reduced.ranks)), c(0, 100), c(1,nrow(reduced.ranks)),c(0,nrow(ranks))),
	xaxis.cex = 0.7,
	xlab.label = c('', '-log10(p)', '', 'Euclidean distance', '', ''),
	ylab.label = c("\t\t\t", "\t\t\t", "\t\t\t", 'Top 5 ranks   Top ranks   '),
	ylab.padding = 1,
	ylab.cex = 0.7,
	bottom.padding = 1,
	xlab.cex = 1,
	yaxis.cex = 0.6,
	x.relation = 'sliced',
	legend = run.legend,
	resolution = 500
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
