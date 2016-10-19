generate.plot.covariates <- function(cov.info) {
	# check input
	avail.covs <- c('type', 'cartridge', 'CodeClass')
	if (!all(names(cov.info) %in% avail.covs)) {
		stop(paste("Sorry, currently only accepting covariates:", paste(avail.covs, collapse = " ")));
		}
	if (class(cov.info) != 'data.frame') {
		stop("cov.info must be of class 'data.frame'!");
		}

	# set up
	covs <- cov.info;
	cols <- list();

	for (i in names(cov.info)) {
		covs[,i] <- as.factor(cov.info[, i]);

		if (i == 'type') {
			cols[[i]] <- colours()[c(507, 532)[1:nlevels(covs[[i]])]];
		} else if (i == 'cartridge') {
			cols[[i]] <- BoutrosLab.plotting.general::colour.gradient('purple', nlevels(covs[[i]]));
		} else if (i == 'CodeClass') {
			cols[[i]] <- BoutrosLab.plotting.general::default.colours(nlevels(covs[[i]]));
			} 
		}
	stopifnot(length(cols) == ncol(covs));

	# create cov rectangles
	output <- list();
	for (i in names(covs)) {
		stopifnot(class(covs[[i]]) == 'factor', nlevels(covs[[i]]) == length(cols[[i]]));
		output[[i]] <- list(
			col = 'transparent',
			fill = cols[[i]][as.numeric(covs[, i])]
			);
		}
	names(output) <- rep('rect', ncol(covs));

	return(output);
	}
