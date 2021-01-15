################################
## Author     : Can Kocagil  ###
################################


rm(list = ls())

data <- read.table("http://www.bixsolutions.net/profiles.csv",sep = ",", header = T)

rangescale <- function(X){
  
  
  Xmax = apply(X,2,max)
  Xscaled <- scale(X,scale = Xmax,center =F)
  
  return(Xscaled)
}



result = prcomp(rangescale(data),center = F)

result$rotation
