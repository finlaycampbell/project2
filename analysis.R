source("R/functions.R")

library(igraph)
library(outbreaker)

temp.outbreak <- simOutbreak(1.5,dgamma(1:10,4,1),n.hosts=20)
plot.simOutbreak(temp.outbreak)

temp.CTD <- simCTD(temp.outbreak,eps=1,avg.contact=4,N=20,plot=TRUE,infec.reported=TRUE)