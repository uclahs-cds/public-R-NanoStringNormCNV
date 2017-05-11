get.ari <- function(data.to.cluster, feature, is.discrete = TRUE){
	if (length(unique(feature)) < 2) {
		flog.warn("..Skipping calculation: number of feature levels must be 2 or more!");
		ari <- NA;
	} else if (is.discrete) {
		# if data is discrete, use Jaccard and Ward.
		dist.matrix <- BoutrosLab.dist.overload::dist(x = t(data.to.cluster), method = 'jaccard');
		corr.hc 	<- hclust(dist.matrix, method = 'ward.D');
		prediction 	<- cutree(corr.hc, length(unique(feature)));
		ari			<- adjustedRandIndex(prediction, factor(feature, levels = unique(feature)));
	} else {
		# if data is continuous, use Pearson and complete clustering.
		corr.matrix <- cor(data.to.cluster, use = 'all.obs', method = 'pearson');
		corr.dist 	<- as.dist(1 - cor(corr.matrix, use = 'pairwise', method = 'pearson'));
		corr.hc 	<- hclust(corr.dist, method = 'complete');
		feature 	<- factor(feature, levels = unique(feature));
		prediction 	<- cutree(corr.hc, length(levels(feature)));
		ari 		<- adjustedRandIndex(prediction, feature);
		}

	return(ari);
	}
