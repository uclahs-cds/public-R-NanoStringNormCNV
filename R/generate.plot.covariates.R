generate.plot.covariates <- function(plotting.data, sample.covariates = NULL, gene.covariates = NULL) {
	# check input
	available.samp.covs <- c('type', 'cartridge');
	available.gene.covs <- c('CodeClass');

	if (!all(names(sample.covariates) %in% c(available.samp.covs, 'SampleID'))) {
		stop(paste("Currently, only accepting sample covariates:", paste(available.samp.covs, collapse = " ")));
		}
	if (!all(names(gene.covariates) %in% c(available.gene.covs, 'Name'))) {
		stop(paste("Currently, only accepting gene covariates:", paste(available.gene.covs, collapse = " ")));
		}
	if (!any(c(class(sample.covariates), class(gene.covariates)) %in% c('matrix', 'data.frame', 'NULL'))) {
		stop("Covariate input must be of class 'data.frame'!");
		}

	# check and set up sample covariates
	if (!is.null(sample.covariates)) {
		# check completeness and order
		if (! all( colnames(plotting.data) %in% sample.covariates$SampleID )) {
			stop("Must provide covariate information for every sample!");
		} else {
			sample.covariates <- sample.covariates[match(
				colnames(plotting.data),
				sample.covariates$SampleID
				),];
			}

		sample.covariates <- sample.covariates[, !(names(sample.covariates) == 'SampleID'), drop = FALSE];
		rownames(sample.covariates) <- NULL;
		}

	# check and set up gene covariates
	if (!is.null(gene.covariates)) {
		# check completeness and order
		if (! all( rownames(plotting.data) %in% gene.covariates$Name )) {
			stop("Must provide covariate information for every gene!");
		} else {
			gene.covariates <- gene.covariates[match(
				rownames(plotting.data),
				gene.covariates$Name
				),];
			}

		gene.covariates <- gene.covariates[, !(names(gene.covariates) == 'Name'), drop = FALSE];
		rownames(gene.covariates) <- NULL;
		}

	# set up
	samp.covs <- sample.covariates;
	gene.covs <- gene.covariates;
	samp.cols <- list();
	gene.cols <- list();

	for (i in names(samp.covs)) {
		samp.covs[,i] <- as.factor(samp.covs[, i]);

		if (i == 'type') {
			samp.cols[[i]] <- colours()[c(507, 532)[1:nlevels(samp.covs[[i]])]];
		} else if (i == 'cartridge') {
			samp.cols[[i]] <- BoutrosLab.plotting.general::colour.gradient('purple', nlevels(samp.covs[[i]]));
			} 
		}
	stopifnot(length(samp.cols) == ncol(samp.covs));

	for (i in names(gene.covs)) {
		gene.covs[,i] <- as.factor(gene.covs[, i]);

		if (i == 'CodeClass') {
			gene.cols[[i]] <- BoutrosLab.plotting.general::default.colours(nlevels(gene.covs[[i]]));
			} 
		}
	stopifnot(length(gene.cols) == ncol(gene.covs));

	# create cov rectangles
	samp.output <- list();
	gene.output <- list();
	
	if (!is.null(samp.covs)) {
		for (i in names(samp.covs)) {
			stopifnot(class(samp.covs[[i]]) == 'factor', nlevels(samp.covs[[i]]) == length(samp.cols[[i]]));
			samp.output[[i]] <- list(col = 'transparent', fill = samp.cols[[i]][as.numeric(samp.covs[,i])]);
			}
		names(samp.output) <- rep('rect', ncol(samp.covs));
		}

	if (!is.null(gene.covs)) {
		for (i in names(gene.covs)) {
			stopifnot(class(gene.covs[[i]]) == 'factor', nlevels(gene.covs[[i]]) == length(gene.cols[[i]]));
			gene.output[[i]] <- list(col = 'transparent', fill = gene.cols[[i]][as.numeric(gene.covs[,i])]);
			}
		names(gene.output) <- rep('rect', ncol(gene.covs));
		}

	return(list(sample = samp.output, gene = gene.output));
	}
