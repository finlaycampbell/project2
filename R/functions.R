#' project2 - Functions
#' 
#' @author F. Campbell
#' 
#' @param R0 the effective reproduction number
#' @param gamma_a the shape value of the gamma function used for w.dens
#' @param gamma_b the rape values of the gamma function used for w.dens

plot.simOutbreak <- function(R0=2,gamma_a=4,gamma_b=1){
  temp.sim = simOutbreak(R0,dgamma(1:10,4,1))
  if(length(temp.sim$id)==1) return("No outbreak")
  temp.net = graph.data.frame(cbind(temp.sim$ances[-1],temp.sim$id[-1]),directed=T)
  return(plot(temp.net))
}

#Returns the genetic pseudo-likelihood using the number of genetic differences, length of comparable genome,
#number of generations and mutation rate
mu.func <- function(d,l,k,u){
  return(u^(d)*(1-u)^(k*l - d))
}

#Plots the likelihood given fixed d, l and k across a range of values for u, returning most likely u value
plot.mu.func <- function(d,l,k){
temp.seq = seq(0,1,0.01)
temp.output = mu.func(d,l,k,temp.seq)
plot(temp.seq,temp.output,type='l',xlab="Mutation rate",ylab="Likelihood")
return(temp.seq[which.max(temp.output)])
}

#Returns the probability of infection from simOutbreak given R0 and n (summing over w(t-ti))
p.inf <- function(R0,n){
  return(1 - exp(1)^(-R0/n))
}