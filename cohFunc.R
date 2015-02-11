cohFunc <- function(x,y) {
  
  m  <- matrix(c(x,y),ncol=2)
  s  <- SDF(m, method="multitaper", n.taper=7)
  sm <- as.matrix(s)
  
  #csd <- Re(sm[,2])
  csn <- abs(sm[,2] * Conj(sm[,2]))
  coh <- csn / Re(sm[,1] * sm[,3])

  return(mean(coh))
}