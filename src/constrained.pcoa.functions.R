require(calibrate)
require(Biostrings)
require(ape)

# function from 
# https://github.com/hzi-bifo/BarleyMeta
# please cite "Davide Bulgarelli, Ruben Garrido-Oter, Philipp C. Münch, Aaron Weiman, Johannes Dröge, Yao Pan, Alice C. McHardy, Paul Schulze-Lefert Structure and Function of the Bacterial Root Microbiota in Wild and Domesticated Barley. Cell Host and Microbe, 2015" if you use this functions

variability_table <- function(cca){
  
  chi <- c(cca$tot.chi,
           cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
  
}

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig/eig_tot
  return(var_propdf)
}

pca_var_props <- function(cca){
  
  eig_tot <- sum(cca$CA$eig)
  var_propdf <- cca$CA$eig/eig_tot
  return(var_propdf)
}

cca_ci <- function(cca, permutations=5000){
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
  return(ci)
  
}


plot_cap <- function(p, d, col_var, pch_var, col_comp, pch_comp,
                     shapes, colors, file_name, svg=F, constraint, 
                     ci, var_tbl, cap_var, perm_anova, cex.main=0.8, cex=0.8) {
  
  col <- rep("black", dim(p)[1])
  pch <- rep(19, dim(p)[1])
  for (i in 1:length(col_var)){
    index <- grep(col_var[i], d[rownames(p), col_comp[i]]) 
    col[index] <- colors[i]
  }
  for (i in 1:length(pch_var)){
    index <- grep(pch_var[i], d[rownames(p), pch_comp[i]]) 
    pch[index] <- shapes[i]
  }
  xlab <- paste("Constrained PCoA 1 (", format(cap_var[1]*100, digits=4), " %)", sep="")
  ylab <- paste("Constrained PCoA 1 (", format(cap_var[2]*100, digits=4), " %)", sep="")
  main <- paste(constraint, ": [", format(var_tbl["constrained", "proportion"]*100, digits=2),
                "% of variance; P < ", format(perm_anova[1,4], digits=2),
                "; 95% CI = ", format(ci[1]*100, digits=2), 
                "%, ", format(ci[2]*100, digits=2), "%]", sep="") 
  if(svg) svg(file=file_name) else pdf(file=file_name)
  plot(p, col=col, pch=pch, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex=cex)
  abline(v=0, h=0, lty="dotted")
  dev.off() 
}