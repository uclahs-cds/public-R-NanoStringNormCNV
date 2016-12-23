get.sample.specific.cna.thresholds <- function(ratios, percent = 0.8) {
	# perform data check
	if(length(ratios) < 2){
		stop("You must specify a data vector of length > 1");
		} 

	# calculate thresholds based on a kernel density method
    dx <- density(x = ratios, kernel = 'gaussian', na.rm = TRUE);
    dn <- cumsum(dx$y) / sum(dx$y);
    li <- which(dn >= (1 - percent) / 2)[1];
    ui <- which(dn >= 1 - (1 - percent) / 2)[1];

    thresholds <- c(
    	lower <- dx$x[li],
    	upper <- dx$x[ui]
    	);

	thresh.df <- data.frame(
		Loss = thresholds[1], 
		Gain = thresholds[2]
		);  

	return(thresh.df);
	}
