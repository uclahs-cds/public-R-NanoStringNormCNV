
get.ari <- function(data, feature, discrete.data = FALSE){
	# if data isn't continuous, use Jaccard and Ward combo
	if(discrete.data){
		dist.matrix <- dist(x=data, method = 'jaccard');
		corr.hc <- hclust(dist.matrix, method = 'ward');
		predicted.feature <- cutree(corr.hc, length(unique(feature)));
		ari <- adjustedRandIndex(predicted.feature, factor(feature, levels = unique(feature)));
		}
	else{
		cor.matrix <- cor(data, use = 'all.obs', method = 'pearson');
		cor.dist <- as.dist(1-cor(cor.matrix, use = 'pairwise', method = 'pearson'));
		cor.dendo <- hclust(cor.dist, method = 'complete'); # 
		feature <- factor(feature, levels = unique(feature));
		predicted.feature <- cutree(cor.dendo, length(levels(feature)));
		ari <- adjustedRandIndex(predicted.feature, feature);
		}

	return(ari);
	}