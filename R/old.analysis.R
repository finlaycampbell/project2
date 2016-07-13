source("R/functions.R")

library(visNetwork)
library(outbreaker)



###################################################################################################
######### A rudimentary means of inferring the transmission tree using by creating random #########
######### trees from simOutbreak and accepting the new tree if it's omega.3 likelihood is #########
######### higher than the current likelihood. Very inefficient and slow but demonstrates  #########
######### that our omega.3 likelihood can inform of the accuracy of our proposal tree     #########
###################################################################################################


#overlay tranmission network

actual.outbreak <- simOutbreak(1.5,dgamma(1:10,4,1),n.hosts=10)
actual.CTD <- simCTD(actual.outbreak,eps=0.7,chi=1,ksi=0.0,plot=TRUE)

sim.ll <- function(temp.eps=0.8){
  prop.outbreak <- simOutbreak(1.5,dgamma(1:10,4,1),n.hosts=5)
  prop.contact <- cbind(prop.outbreak$ances,prop.outbreak$id)
  import <- is.na(prop.contact[,1])
  prop.contact <- prop.contact[-import,]
  if(is.null(nrow(prop.contact))) return(list(outbreak=prop.outbreak,ll=0)) else if(nrow(prop.contact) == 0) return(list(outbreak=prop.outbreak,ll=prop.ll))
  prop.ll <- ll.CTD(actual.CTD,prop.contact,temp.eps)
  return(list(outbreak=prop.outbreak,ll=prop.ll))
}

current.it = sim.ll()

for(i in 1:iterations){
  temp.it = sim.ll()
  if(temp.it$outbreak$n<0.9*actual.outbreak$n) next
  if(temp.it$ll > current.it$ll) current.it <- temp.it
}

simCTD(current.it$outbreak,plot=TRUE)