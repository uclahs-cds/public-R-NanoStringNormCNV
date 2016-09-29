
get.ari <- function(data, feature, discrete.data = FALSE){
	# if data isn't continuous, use Jaccard and Ward in combination
	if (discrete.data) {
		dist.matrix <- dist(x = data, method = 'jaccard');
		corr.hc 	<- hclust(dist.matrix, method = 'ward.D');
		prediction 	<- cutree(corr.hc, length(unique(feature)));
		ari			<- adjustedRandIndex(prediction, factor(feature, levels = unique(feature)));
	} else {
		corr.matrix <- cor(data, use = 'all.obs', method = 'pearson');
		corr.dist 	<- as.dist(1 - cor(corr.matrix, use = 'pairwise', method = 'pearson'));
		corr.dendo 	<- hclust(corr.dist, method = 'complete');
		feature 	<- factor(feature, levels = unique(feature));
		prediction 	<- cutree(corr.dendo, length(levels(feature)));
		ari 		<- adjustedRandIndex(prediction, feature);
		}

	return(ari);
	}
