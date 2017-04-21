pairwisePermanova <- function(dist, map){
  require(vegan)
  require(pander)
  F.Model <- NULL; R2 <- NULL; pairs <- NULL;  p.value <- NULL
  comp.list = combn(unique(map),2) # list of all possible one vs. on combinations 
  for(comp.index in 1:ncol(comp.list)){
    map.index <- which(map == comp.list[1,comp.index] | map == comp.list[2,comp.index])
    subset.dist <- dist[map.index, map.index]  # distance matrix reduced to samples that are tested
    subset.map <- map[map.index] # mapping file reduced to samples that are tested
    model <- adonis(as.dist(subset.dist) ~ subset.map, permutations = 5000) # fit the model
    pairs <- c(pairs,paste(comp.list[1,comp.index],'vs',comp.list[2,comp.index])) 
    F.Model <- c(F.Model,model$aov.tab[1,4])
    R2 <- c(R2,model$aov.tab[1,5])
    p.value <- c(p.value,model$aov.tab[1,6])
  }
  p.adjusted <- p.adjust(p.value,method="fdr")
  p.stars <- add.significance.stars(p.adjusted)
  pairw.res <- data.frame(pairs, F.Model, R2, p.value, p.adjusted, p.stars)
  return(pairw.res)
}

pairwiseCap<- function(dist, map){
  require(vegan)
  require(pander)
  F.Model <- NULL; R2 <- NULL; pairs <- NULL;  p.value <- NULL
  comp.list = combn(unique(map),2) # list of all possible one vs. on combinations 
  for(comp.index in 1:ncol(comp.list)){
    map.index <- which(map == comp.list[1,comp.index] | map == comp.list[2,comp.index])
    subset.dist <- dist[map.index, map.index]  # distance matrix reduced to samples that are tested
    subset.map <- map[map.index] # mapping file reduced to samples that are tested
    cap <- capscale(as.dist(subset.dist) ~ subset.map) # fit the model
    permanova <- anova.cca(cap, perm = 5000)
    pairs <- c(pairs,paste(comp.list[1,comp.index],'vs',comp.list[2,comp.index])) 
    p.value <- c(p.value, permanova$`Pr(>F)`[1])
  }
  p.adjusted <- p.adjust(p.value,method="fdr")
  p.stars <- add.significance.stars(p.adjusted)
  pairw.res <- data.frame(pairs, p.value, p.adjusted, p.stars)
  return(pairw.res)
}