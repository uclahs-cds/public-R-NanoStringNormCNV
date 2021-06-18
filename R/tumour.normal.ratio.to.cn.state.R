tumour.normal.ratio.to.cn.state <- function(ratios, thresholds) {
	# helper function to call five-state data:
	# -2, -1, 0, 1, 2; hom del, het del, neutral, single copy gain, 2+ copy gain
	call.5state <- function (thresholds, ratio) {
		if (is.na(ratio)) {
			return(NA);
		} else if (ratio <= thresholds[1]) {
			return(-2);
		} else if (ratio <= thresholds[2]) {
			return(-1);
		} else if (ratio >= thresholds[4]) {
			return(2);
		} else if (ratio >= thresholds[3]) {
			return(1);
		} else {
			return(0);
			}
		}

	# iterate through each segment
	cn.state <- rep(0, length(ratios));
	for (seg in 1:length(ratios)) {
		cn.state[seg] <- call.5state(thresholds, ratios[seg]);
		}

	return(cn.state);
	}
