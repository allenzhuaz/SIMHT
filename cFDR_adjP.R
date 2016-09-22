#' @param pval The structural p-values, the type is a list.  
#' @param t The thresholds determine the families are selected or not, also affects conditional p-value within families
#' @param method p-value combining methods

cFDR.cp.adjust <- function(pval, t, method=c("Fisher", "Stouffer", "minP")){
  
  method <- match.arg(method)
  if (length(t)<length(pval)){t <- rep(t,len=length(pval))}
  condPval <- vector("list", length(pval))
  adjcondP <- vector("list", length(pval))
 
  if (method=="Fisher"){
    for (i in which(lapply(pval,prod)<=exp(-t/2)) ){
      condPval[[i]] <- integer(length(pval[[i]]))
      for (j in 1:length(pval[[i]])){
        if (length(pval[[i]])==1){condPval[[i]][j] <- pval[[i]][j]/exp(-t[i]/2)}
        else {
          condPval[[i]][j] <- pval[[i]][j]/min(c(exp(-t[i]/2)/prod(pval[[i]][-j]),1))
        }
      }
      adjcondP[[i]] <- p.adjust(condPval[[i]], method = "BH")
    }
  }
  
  if (method=="Stouffer"){
    stfun <- function(x){sum(qnorm(1-x))/sqrt(length(x))}
    for (i in which(sapply(pval,stfun)>=t) ){
      condPval[[i]] <- integer(length(pval[[i]]))
      for (j in 1:length(pval[[i]])){
        if (length(pval[[i]])==1){condPval[[i]][j] <- pval[[i]][j]/(1-pnorm(t[i]))}
        else {
        condPval[[i]][j] <- pval[[i]][j]/(1-pnorm(sqrt(length(pval[[i]]))*t[i]-sum(qnorm(1-pval[[i]][-j]))))
        }
      }
      adjcondP[[i]] <- p.adjust(condPval[[i]], method = "BH")
    }
  }
  
  else if (method=="minP"){
    for (i in which(sapply(pval,min)<=t) ){
      condPval[[i]] <- integer(length(pval[[i]]))
      for (j in 1:length(pval[[i]])){
        if (length(pval[[i]])==1){condPval[[i]][j] <- pval[[i]][j]/t[i]}
        else {
          condPval[[i]][j] <- pval[[i]][j]/ifelse(min(pval[[i]][-j]) > t[i], t[i] ,1)
        }
      }
      adjcondP[[i]] <- p.adjust(condPval[[i]], method = "BH")
    }
  }
    return(adjcondP)
}