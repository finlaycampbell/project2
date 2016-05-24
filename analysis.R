source("R/functions.R")

library(igraph)
library(visNetwork)
library(outbreaker)



temp.outbreak <- simOutbreak(1.5,dgamma(1:10,4,1),n.hosts=20)
plot.simOutbreak(temp.outbreak)

temp.CTD <- simCTD(temp.outbreak,eps=1,chi=1,ksi=0,plot=TRUE)


actual.outbreak <- simOutbreak(1.5,dgamma(1:10,4,1),n.hosts=10)
plot.simOutbreak(actual.outbreak)
actual.CTD <- simCTD(actual.outbreak,eps=1,avg.contact=3,N=300,plot=TRUE)

iterations = 50

for(i in 1:iterations){
  prop.outbreak <- simOutbreak(3,dgamma(1:10,4,1),n.hosts=10)
  prop.contact <- cbind(prop.outbreak$ances[-1],prop.outbreak$id[-1])
  prop.contact = prop.contact[-which(is.na(prop.contact[,1])),]
  if(nrow(prop.contact)==0) next
  prop.ll = ll.CTD(actual.CTD,prop.contact,0.8)
  print(prop.ll)
}