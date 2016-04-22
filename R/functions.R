#DEFINING FUNCTIONS

##Visualising simOutbreak
plot.simOutbreak <- function(R0=2,gamma_a=4,gamma_b=1){
  library(igraph)
  library(network)
  
  temp.sim = simOutbreak(1.8,dgamma(1:10,4,1))
  temp.net = graph.data.frame(cbind(temp.sim$ances,temp.sim$id),directed=T)
  return(plot(temp.net))
}

##Returns the genetic pseudo-likelihood using the number of genetic differences, length of comparable genome,
##number of generations and mutation rate
mu.func = function(d,l,k,u){
  return(u^(d)*(1-u)^(k*l - d))
}

##Plots the likelihood given fixed d, l and k across a range of values for u, returning most likely u value
plot.mu.func = function(d,l,k){
temp.seq = seq(0,1,0.01)
temp.output = mu.func(d,l,k,temp_seq)
plot(temp.seq,temp/output,type='l')
return(temp.seq[which.max(temp.output)])
}