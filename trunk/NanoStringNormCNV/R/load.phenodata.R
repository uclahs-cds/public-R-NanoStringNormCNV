load.phenodata <- function(fname, separator = 'comma') {
	# read in data (by file type)
	if (separator == "comma") {
		phenodata <- read.csv(fname, stringsAsFactors = FALSE);
	} else if (separator == "tab") {
		phenodata <- read.delim(fname, stringsAsFactors = FALSE);
		}

	# check for mandatory column names
	required.cols <- c(
		"SampleID",
		"Patient",
		"Name",
		"Cartridge",
		"Type",
		"ReferenceID",
		"HasReplicate"
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
	if (any(duplicated(phenodata$SampleID))) {
		stop("Sample IDs must be unique!");
		}

	# check that sample Names match across replicates
	unmatched.reps <- c();
	phenodata.for.reps <- phenodata[phenodata$HasReplicate == 1,];
	for (i in unique(phenodata.for.reps$Name)) {
		if (length(which(phenodata.for.reps$Name == i)) < 2) {
			unmatched.reps <- c(unmatched.reps, paste0("\n\t", i));
			}
		}

	if (length(unmatched.reps) > 0) {
		flog.warn(paste0(
			"Sample replicates must have matching sample names (see column 'Name')! ",
			"Cannot identify replicate(s) for the following samples:",
			paste0(unmatched.reps, collapse = "")
			));
		}
	
	# check Cartridge values are numeric
	if (!is.numeric(phenodata$Cartridge)) {
		flog.warn("Cartridge values must be numeric");
		phenodata$Cartridge <- as.numeric(phenodata$Cartridge);
		if (any(is.na(phenodata$Cartridge))) { stop("Unable to convert cartridge values to numeric!") }
		}

	# check Type values
	phenodata$Type[tolower(phenodata$Type) == 'tumour' | tolower(phenodata$Type) == 'tumor'] <- 'Tumour';
	phenodata$Type[tolower(phenodata$Type) == 'reference'] <- 'Reference';
	if (!all(phenodata$Type == 'Tumour' | phenodata$Type == 'Reference')) {
		stop("Column 'Type' must contain only 'Tumour' or 'Reference'!");
		}

	# check reference sample information
	ref.tumour <- phenodata[phenodata$Type == 'Tumour',]$ReferenceID;
	ref.normal <- phenodata[phenodata$Type == 'Reference',]$ReferenceID;

	for (i in which(!is.na(phenodata$ReferenceID) & phenodata$ReferenceID != 'missing')) {
		if (!(phenodata$ReferenceID[i] %in% phenodata$SampleID)) {
			stop(paste0(
				"Cannot identify reference sample ", phenodata[i,]$ReferenceID,
				" for tumour sample ", phenodata[i,]$SampleID, "!"
				 ));
			}
		}

	if (length(which(phenodata$Type == 'Reference')) < 1) {
		flog.warn("Column 'Type' contains no reference samples: will be unable to call CNAs downstream!");
		}

	if (any(!(ref.tumour[ref.tumour != 'missing'] %in% phenodata$SampleID))) {
		stop(paste0(
			"Column 'ReferenceID' must contain sample IDs of matched normal samples. ",
			"For tumour sample without matched normals, put 'missing'."
			));
		}

	if (any(!is.na(ref.normal))) {
		stop(paste0(
			"Column 'ReferenceID' must contain sample IDs of matched normal samples. ",
			"For normal samples, put NA."
			));
		}

	# check replicate information
	if (!all(phenodata$HasReplicate == 0 | phenodata$HasReplicate == 1)) {
		stop(paste0(
			"Column 'HasReplicate' must contain only the following values:",
			"\n\t1 for sample with replicate",
			"\n\t0 for sample without replicate"
			));
		}

	# check for sex information
	if ("sex" %in% tolower(names(phenodata))) {
		names(phenodata)[which("sex" == tolower(names(phenodata)))] <- "Sex";

		phenodata$Sex[tolower(phenodata$Sex) %in% c("female", "f")] <- "F";
		phenodata$Sex[tolower(phenodata$Sex) %in% c("male", "m")]   <- "M";

		if (any(na.omit(phenodata$Sex) != "M" & na.omit(phenodata$Sex) != "F")) {
			stop("Column 'Sex' must contain only the following values: M, F, NA");
			}

		# all samples belonging to one patient should have the same sex!
		for (i in unique(phenodata$Patient)) {
			if (length(unique(phenodata[phenodata$Patient == i,]$Sex)) > 1) {
				stop(paste0("More than one sex assigned to patient ", i, "!"));
				}
			}

	} else {
		flog.warn("Column 'Sex' not provided: unable to call CNAs on sex chromosomes downstream!");
		phenodata$Sex <- NA;
		}

	# ensure 'SampleID' and 'Name' values don't contain special characters
	pattern <- "[^a-zA-Z0-9\\.]";

	if (any(grepl(pattern, phenodata$SampleID))) {
		flog.info("'SampleID' values should only contain alphanumeric characters. Replacing non-alphanumeric with '.'");
		phenodata$ReferenceID[phenodata$ReferenceID %in% phenodata$SampleID] <- gsub(pattern, "\\.", phenodata$ReferenceID[phenodata$ReferenceID %in% phenodata$SampleID]);
		phenodata$SampleID <- gsub(pattern, "\\.", phenodata$SampleID);
		}

	if (any(grepl(pattern, phenodata$Name))) {
		flog.info("'Name' values should only contain alphanumeric characters. Replacing non-alphanumeric with '.'");
		phenodata$Name <- gsub(pattern, "\\.", phenodata$Name);
		}

	return(phenodata);
	}
