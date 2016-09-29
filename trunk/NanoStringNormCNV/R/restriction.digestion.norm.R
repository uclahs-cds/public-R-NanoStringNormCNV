
### Step 2: Restriction fragmentation normalization
restriction.fragmentation.norm <- function(nano.df){
	#--- restriction fragmentation controls ----------------------------------------------------------#
	# restriction controls evaluation
	nano.restr <- nano.df[nano.df$CodeClass == 'RestrictionSite',];
	nano.restr <- nano.restr[with(nano.restr, order(Name)),];

	# remove CodeClass, Name, Accession from colnames
	grep.cols <- colnames(nano.restr)[!colnames(nano.restr) %in% c('CodeClass', 'Name', 'Accession')];

	# create a new data-frame
	nano.restr.avg <- as.data.frame(
		matrix(
			ncol = length(grep.cols),
			nrow = 2,
			dimnames = list(
				c('A+B', 'C+D'),
				grep.cols
				)
			),
		stringsAsFactors = FALSE
		);

	# grep A and B names
	grep.AB <- c(
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE+A', fixed = TRUE),
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE+B', fixed = TRUE)
		);

	grep.CD <- c(
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE-C', fixed = TRUE),
		grep(x = nano.restr$Name, pattern = 'RESTRICTIONSITE-D', fixed = TRUE)
		);

	# calculate the count average for A,B and C,D
	nano.restr.avg['A+B',] <- apply(X = nano.restr[grep.AB,grep.cols], MARGIN = 2, FUN = mean);
	nano.restr.avg['C+D',] <- apply(X = nano.restr[grep.CD,grep.cols], MARGIN = 2, FUN = mean);

	# if low counts is observed in RESTRICTIONSITE-C and RESTRICTIONSITE-D (<200), throw a warning
	if (!all(nano.restr.avg['C+D',] >= 200)) {

		# list out the sample names
		flog.warn('Low count number!');
		low.count <- paste(colnames(nano.restr.avg)[nano.restr.avg['C+D',] < 200], sep = '\n');
		print(low.count);
		}

	# check if average of A and B and average of C and D has a 10-fold difference in counts
	which.low <- NULL;
	if (!all((nano.restr.avg['C+D',]/nano.restr.avg['A+B',]) > 55)) {

		# gives a warning, record the ratio and see if these should be removed
		which.low <- which((nano.restr.avg['C+D',]/nano.restr.avg['A+B',]) <= 55);

        # get the data-frame
        low.ratio <- nano.restr.avg[, which.low, drop = FALSE];

        # add a new row with the ratio
        low.ratio['ratio',] <- low.ratio['C+D',]/low.ratio['A+B',];

        # throw out a warning
        flog.warn('Low Fold-Change!', capture = TRUE, low.ratio);
        }

		# plot the ratios
		NanoStringNormCNV::make.restr.digest.plot(nano.restr.avg, which.low);

		out.df <- data.frame(
			CD = as.numeric(t(nano.restr.avg['C+D',])),
			AB = as.numeric(t(nano.restr.avg['A+B', ])),
			ratio = as.numeric(t(nano.restr.avg['C+D',]/nano.restr.avg['A+B',]))
			);
		rownames(out.df) <- colnames(nano.restr.avg);
		return(out.df);
	}
