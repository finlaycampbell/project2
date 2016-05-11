source("R/functions.R")

library(ggplot2)
library(igraph)
library(network)
library(outbreaker)

##Constructing simulated genomes and collection dates for use by outbreaker
test.dna = matrix(c("C","C","C","C","C","C","C","C","C",
                    "C","G","C","C","C","C","C","C","C",
                    "C","G","C","T","C","C","C","C","C",
                    "C","G","C","T","C","A","C","C","C",
                    "C","G","C","T","C","A","T","C","C"),nrow=5)

test.date = c(0,4,7,10,14)
test.w.dens = dgamma(1:11,2,0.5)
test.output = outbreaker(test.dna,test.date,w.dens=test.w.dens)

temp.sim <- simOutbreak(1.5,dgamma(1:10,4,1),n.hosts=100)
report.contact <- simContact(temp.sim,eps=0.5,avg.contact=8,N=1000)
plot(graph.data.frame(report.contact))

true.contact.reported(temp.sim,report.contact)