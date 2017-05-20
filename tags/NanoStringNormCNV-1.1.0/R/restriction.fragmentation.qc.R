restriction.fragmentation.qc <- function(raw.data){
	# only keep the rows containing restriction fragmentation controls
	nano.restr <- raw.data[raw.data$CodeClass == 'RestrictionSite',];
	nano.restr <- nano.restr[with(nano.restr, order(Name)),];

	grep.cols <- colnames(nano.restr)[!colnames(nano.restr) %in% c('CodeClass', 'Name', 'Accession')];

	nano.restr.avg <- as.data.frame(
		matrix(
			ncol = length(grep.cols),
			nrow = 2,
			dimnames = list(c('A+B', 'C+D'), grep.cols)
			),
		stringsAsFactors = FALSE
		);

	grep.AB <- c(
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE+A', fixed = TRUE),
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE+B', fixed = TRUE)
		);

	grep.CD <- c(
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE-C', fixed = TRUE),
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE-D', fixed = TRUE)
		);

	# calculate mean counts probes for A+B and C+D
	nano.restr.avg['A+B',] <- apply(X = nano.restr[grep.AB,grep.cols], MARGIN = 2, FUN = mean);
	nano.restr.avg['C+D',] <- apply(X = nano.restr[grep.CD,grep.cols], MARGIN = 2, FUN = mean);

	# check if low counts (<200) are observed in probes C+D
	if (!all(nano.restr.avg['C+D',] >= 200)) {

		# list out the sample names
		flog.warn('Low count number!');
		low.count <- paste(colnames(nano.restr.avg)[nano.restr.avg['C+D',] < 200], sep = '\n');
		print(low.count);
		}

	# check if there is a 10-fold difference between A+B probe mean counts and C+D probe mean counts
	which.low <- NULL;
	if (!all((nano.restr.avg['C+D',]/nano.restr.avg['A+B',]) > 10)) {

		which.low <- which((nano.restr.avg['C+D',]/nano.restr.avg['A+B',]) <= 10);
        low.ratio <- nano.restr.avg[, which.low, drop = FALSE];
        low.ratio['ratio',] <- low.ratio['C+D',]/low.ratio['A+B',];
        flog.warn('Low Fold-Change!', capture = TRUE, low.ratio);

        }

	# plot the ratios
	NanoStringNormCNV:::make.restriction.fragmentation.plot(nano.restr.avg, which.low);

	# set up output
	out.df <- data.frame(
		CD = as.numeric(t(nano.restr.avg['C+D',])),
		AB = as.numeric(t(nano.restr.avg['A+B', ])),
		ratio = as.numeric(t(nano.restr.avg['C+D',]/nano.restr.avg['A+B',]))
		);
	rownames(out.df) <- colnames(nano.restr.avg);

	return(out.df);
	}
