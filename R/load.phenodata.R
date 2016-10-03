# do something with outliers here?
# check if 'lane', 'tissue' info is used
load.phenodata <- function(fname){
	# read data
	phenodata <- read.csv(fname, stringsAsFactors = FALSE);

	# check column names
	required.cols <- c(
		"SampleID",
		"Patient",
		"Name",
		"cartridge",
		"type",
		"ref.name",
		"has.repl"
		);

	for (i in required.cols) {
		if (!(i %in% names(phenodata))) {
			stop(paste0(
				"Missing annotation column \'", i, "\'.",
				" Please see documentation for more information.\n"
				));
			}
		}

	# check Sample IDs are unique
	if (any(duplicated(phenodata$SampleID)) | any(duplicated(phenodata$Name))) {
		stop("Sample IDs and names must be unique!");
		}

	# check cartridge values are numeric
	if (!is.numeric(phenodata$cartridge)) {
		flog.warn("Cartridge values must be numeric");
		phenodata$cartridge <- as.numeric(phenodata$cartridge);
		if (any(is.na(phenodata$cartridge))) { stop("Unable to convert cartridge values to numeric!") }
		}

	# check type values
	phenodata$type[tolower(phenodata$type) == 'tumour']    <- 'tumour';
	phenodata$type[tolower(phenodata$type) == 'reference'] <- 'reference';
	if (!all(phenodata$type == 'tumour' | phenodata$type == 'reference')) {
		stop("Column 'type' must contain only 'tumour' or 'reference'!");
		}	

	# check reference sample information
	ref.tumour <- phenodata[phenodata$type == 'tumour',]$ref.name;
	ref.normal <- phenodata[phenodata$type == 'reference',]$ref.name;

	if (any(!(ref.tumour[ref.tumour != 'missing'] %in% phenodata$SampleID))) {
		stop(paste0(
			"Column 'ref.name' must contain Sample IDs of matched normal samples. ",
			"For tumour sample without matched normals, put 'missing'."
			));
		}

	if (any(!is.na(ref.normal))) {
		stop(paste0(
			"Column 'ref.name' must contain Sample IDs of matched normal samples. ",
			"For normal samples, put NA."
			));
		}

	# check replicate information
	if (!all(phenodata$has.repl == 0 | phenodata$has.repl == 1)) {
		stop(paste0(
			"Column 'has.repl' must contain the following values:",
			"\n\t1 for sample with replicate",
			"\n\t0 for sample without replicates"
			));
		}

	# ensure 'SampleID' and 'Name' values don't contain special characters
	pattern <- "[^a-zA-Z0-9\\.]";

	if (any(grepl(pattern, phenodata$SampleID))) {
		flog.warn("'SampleID' values should only contain alphanumeric characters. Replacing non-alphanumeric with '.'");
		phenodata$ref.name[phenodata$ref.name %in% phenodata$SampleID] <- gsub(pattern, "\\.", phenodata$ref.name[phenodata$ref.name %in% phenodata$SampleID]);
		phenodata$SampleID <- gsub(pattern, "\\.", phenodata$SampleID);
		}

	if (any(grepl(pattern, phenodata$Name))) {
		flog.warn("'Name' values should only contain alphanumeric characters. Replacing non-alphanumeric with '.'");
		phenodata$Name <- gsub(pattern, "\\.", phenodata$Name);
		}

	# order data
	phenodata <- phenodata[order(phenodata$SampleID),];

	return(phenodata);
	}
