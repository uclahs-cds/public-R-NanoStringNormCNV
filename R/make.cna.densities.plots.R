
make.gene.densities.plot <- function(cnas, fname.stem, xlab){

	### per gene
	gene.list <- lapply(seq(1:nrow(cnas)), function(f) cnas[f,]);
	names(gene.list) <- paste0('gene', 1:length(gene.list));

	create.densityplot(
		x = gene.list,
		filename = generate.filename(fname.stem, 'gene_density', 'png'),
		xlab.label = xlab,
		xlab.cex = 2,
		main = 'Density per gene',
		main.cex = 2,
		add.rectangle = TRUE,
		xleft.rectangle = c(-5, 0.5, 1.5, 2.5, 3.5),
		ybottom.rectangle = c(0,0,0,0,0),
		xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 15),
		ytop.rectangle = c(50, 50, 50, 50, 50),
		col.rectangle = c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick'),
		alpha.rectangle = 0.6,
		lwd = 0.5,
		lty = 'dotted',
		resolution = 600
		);

	### per sample
	sample.list <- lapply(seq(1:ncol(cnas)), function(f) cnas[, f]);
	names(sample.list) <- paste0('sample', 1:length(sample.list));

	create.densityplot(
		x = sample.list,
		filename = generate.filename(fname.stem, 'sample_density', 'png'),
		xlab.label = xlab,
		xlab.cex = 2,
		main = 'Density per sample',
		main.cex = 2,
		add.rectangle = TRUE,
		xleft.rectangle = c(-5, 0.5, 1.5, 2.5, 3.5),
		ybottom.rectangle = c(0,0,0,0,0),
		xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 15),
		ytop.rectangle = c(50, 50, 50, 50, 50),
		col.rectangle = c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick'),
		alpha.rectangle = 0.6,
		lwd = 0.5,
		lty = 'dotted',
		resolution = 600
		);

	### per sample logged
	sample.list <- lapply(seq(1:ncol(cnas)), function(f) log2(cnas[, f]));
	names(sample.list) <- paste0('sample', 1:length(sample.list));

	create.densityplot(
		x = sample.list,
		filename = generate.filename(fname.stem, 'sample_density-logged', 'png'),
		xlab.label = xlab,
		xlab.cex = 2,
		main = 'Density per sample',
		main.cex = 2,
		add.rectangle = TRUE,
		xleft.rectangle = c(-5, 0.5, 1.5, 2.5, 3.5),
		ybottom.rectangle = c(0,0,0,0,0),
		xright.rectangle = c(0.5, 1.5, 2.5, 3.5, 15),
		ytop.rectangle = c(50, 50, 50, 50, 50),
		col.rectangle = c('dodgerblue4', 'lightblue3', 'white', 'lightcoral', 'firebrick'),
		alpha.rectangle = 0.6,
		lwd = 0.5,
		lty = 'dotted',
		resolution = 600
		);
	
	}
