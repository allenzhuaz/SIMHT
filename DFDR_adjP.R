#' Double FDR (DFDR)
#' @param pval The p-value list
#' @param t The threshold for the family level 


DFDR.p.adjust <- function(pval, t){
  condPval <- vector("list", length(pval))
  adjcondP <- vector("list", length(pval))
  for (i in which(p.adjust(sapply(pval,min),method="BH")<=t) ){
    adjcondP[[i]] <- p.adjust(pval[[i]], method = "BH")
  }
  return(adjcondP)
}