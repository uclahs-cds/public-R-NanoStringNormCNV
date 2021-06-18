generate.plot.legend <- function(cov.info) {
	# build legend
	sample.legend <- list();
	for (i in names(cov.info)) {
		colours <- NULL;
		cov.info[[i]] <- as.factor(cov.info[[i]]);

		if (i == 'CodeClass') {
			colours <- BoutrosLab.plotting.general::default.colours(nlevels(cov.info[[i]]));
		} else if (i == 'Type') {
			colours <- colours()[c(507, 532)[1:nlevels(cov.info[[i]])]];
		} else if (i == 'Cartridge') {
			colours <- BoutrosLab.plotting.general::colour.gradient('purple', nlevels(cov.info[[i]]));
		} else {
			stop("can only accept the following covariates: 'Type', 'Cartridge', 'CodeClass'");
			}

		sample.legend[[i]] <- list(
			colours = colours,
			labels = levels(cov.info[[i]]),
			title = bquote(underline(.(gsub("^([[:alpha:]])", "\\U\\1", i, perl = TRUE))))
			);
		}

	names(sample.legend) <- rep('legend', length(sample.legend));

	return(sample.legend);
	}
