
### normalize all data in one batch
normalize.global <- function(raw.data, cc, bc, sc, do.nsn, do.rcc.inv, covs, plot.types = 'all'){
	if(do.nsn){
		nano.norm <- NanoStringNorm(
			x = raw.data[, -c(1:3)],
			CodeCount = cc,
			Background = bc,
			SampleContent = sc,
			round.values = FALSE,
			take.log = FALSE,
			traits = covs,
			anno = raw.data[,1:3]
			);
		normalized.data <- nano.norm$normalized.data;
		colnames(normalized.data)[1] <- 'CodeClass';
		pdf('NanoStringNorm_plots_all.pdf');
		Plot.NanoStringNorm(x=nano.norm, label.best.guess=T, plot.type=qw('cv mean.sd batch.effects norm.factors missing RNA.estimages positive.controls'));
		dev.off();
		}
	else{
		normalized.data <- raw.data;
		}

	if(do.rcc.inv){
		normalized.data <- NanoStringNormCNV::invariant.probe.norm(normalized.data);
		}
	return(normalized.data);
	}


