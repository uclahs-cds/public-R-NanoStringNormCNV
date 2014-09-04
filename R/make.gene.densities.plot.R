
make.gene.densities.plot <- function(cnas, fname.stem, xlab){
	gene.list <- lapply(seq(1:nrow(cnas)), function(f) cnas[f,]);
	names(gene.list) <- paste0('gene', 1:length(gene.list));

	create.densityplot(
		x = gene.list,
		filename = generate.filename(fname.stem, 'gene_density', 'png'),
		xlab.label = xlab,
		xlab.cex = 2,
		main = 'Density per gene',
		main.cex = 2,
		lwd = 0.5,
		lty = 'dotted',
		resolution = 600
		);
	}
