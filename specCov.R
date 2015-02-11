specCov <- function(x,y) {
  
  specX=spectrum(x,plot=F)$spec
  specY=spectrum(y,plot=F)$spec
  
  return(cov(specX,specY))
}