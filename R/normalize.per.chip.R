
### normalize data for each chip separately
normalize.per.chip <- function(pheno, raw.data, cc, bc, sc, oth, do.nsn, do.rcc.inv, covs, plot.types = 'all'){
	nano.parts <- list();
	cartridges <- unique(pheno$cartridge);
	for(chip in 1:length(unique(pheno$cartridge))){
		cur.samples <- which(pheno$cartridge == cartridges[chip]);
		if(do.nsn){
			if(length(unique(covs[cur.samples , 'Type'])) > 1){
				use.covs <- covs[cur.samples , 'Type', drop = FALSE];
			}else{
				use.covs  <- NA;
				}
			nano.parts[[chip]] <- NanoStringNorm(
				x = raw.data[,(cur.samples+3)],
				CodeCount = cc,
				Background = bc,
				SampleContent = sc,
				Other = oth,
				round.values = FALSE,
				take.log = FALSE,
				traits = use.covs,
				anno = raw.data[, 1:3]
				);
			pdf(paste0('NanoStringNorm_plots_chip', chip, '.pdf'));
			Plot.NanoStringNorm(x=nano.parts[[chip]], label.best.guess=T, plot.type=unlist(strsplit("cv mean.sd norm.factors missing RNA.estimates positive.controls", "\\s")));
			dev.off();
			nano.parts[[chip]] <- nano.parts[[chip]]$normalized.data;
			colnames(nano.parts[[chip]])[1] <- 'CodeClass';
		}else{
			nano.parts[[chip]] <- raw.data[,c(1:3, (cur.samples+3))];
			}
		if(do.rcc.inv){
			nano.parts[[chip]] <- NanoStringNormCNV::invariant.probe.norm(nano.parts[[chip]]);
			}
		}
	normalized.data <- cbind(raw.data[, 1:3], do.call(cbind, lapply(nano.parts, function(f) f[,-c(1:3)])));
	return(normalized.data);
	}

