### tumour.normal.ratio.to.cn.state.R ##########################################################
# Copied (and modified) from BoutrosLab.utilities.copynumber::call.cna.states.R.

tumour.normal.ratio.to.cn.state <- function(ratio.data, thresholds) {
	# helper function to call five-state data:
	# 0, 1, 2, 3, 4; hom del, het del, neutral, single copy gain, 2+ copy gain
	call.5state <- function (thresholds, ratio) {
		if (is.na(ratio)) {
			return(NA);
		} else if (ratio <= thresholds[1]) {
			return(0);
		} else if (ratio <= thresholds[2]) {
			return(1);
		} else if (ratio >= thresholds[4]) {
			return(4);
		} else if (ratio >= thresholds[3]) {
			return(3);
		} else {
			return(2);
			}
		}

	# iterate through each segment
	cn.state <- rep(0, nrow(ratio.data));

	for (seg in 1:length(ratio.data$ratio)) {
		cn.state[seg] <- call.5state(thresholds, ratio.data$ratio[seg]);
		}

	return(cn.state);
	}