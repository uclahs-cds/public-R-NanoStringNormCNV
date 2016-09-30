### call.cna.states.R ##############################################################################
# NOTE:
# Each entry in data.list requires columns named 'chr' and 'log2ratio'
# Copied (and modified) from BoutrosLab.utilities.copynumber::call.cna.states.R.

call.cna.states <- function(cna.data, thresholds) {
	# helper function to call five-state data:
	# 0, 1, 2, 3, 4; hom del, het del, neutral, single copy gain, 2+ copy gain
	call.4state <- function (thresholds, log2ratio) {
		if (is.na(log2ratio)) {
			return(NA);
		} else if (log2ratio <= thresholds[1]) {
			return(0);
		} else if (log2ratio <= thresholds[2]) {
			return(1);
		} else if (log2ratio >= thresholds[4]) {
			return(4);
		} else if (log2ratio >= thresholds[3]) {
			return(3);
		} else {
			return(2);
			}
		}

	# iterate through each segment
	cn.state <- rep(0, nrow(cna.data));

	for (seg in 1:length(cna.data$log2ratio)) {
		cn.state[seg] <- call.4state(thresholds, cna.data$log2ratio[seg]);
		}

	return(cn.state);
	}


