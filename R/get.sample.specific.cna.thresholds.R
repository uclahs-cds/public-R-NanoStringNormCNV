### get.sample.specific.cna.thresholds.R ###########################################################
# Function to determine the thresholds to use for calling CNA states.
# Thresholds are based on the standard deviations of the central 50% of probes (up to user to pass
# the correct median/sd).
# Copied (and modified) from BoutrosLab.utilities.copynumber::get.sample.specific.cna.thresholds.R.

get.sample.specific.cna.thresholds <- function(cna.data, percent = 0.8) {

	# perform data check
	if(length(cna.data) < 2){
		stop("You must specify a data vector of length > 1");
		} 

	# calculate thresholds based on a kernel density method
    dx <- density(x = cna.data, kernel = 'gaussian', na.rm = TRUE);
    dn <- cumsum(dx$y) / sum(dx$y, na.rm = TRUE);
    li <- which(dn >= (1 - percent) / 2)[1];
    ui <- which(dn >= 1 - (1 - percent) / 2)[1];

    thresholds <- c(
    	lower <- dx$x[li],
    	upper <- dx$x[ui]
    	);

	thresh.mat <- data.frame(
		Loss   = thresholds[1], 
		Gain   = thresholds[2]
		);  

	return(thresh.mat);
	}
