#' project2 - Functions
#' 
#' 
#' 
#' @author F. Campbell
#' 
#' @param R0 the effective reproduction number
#' @param gamma_a the shape value of the gamma function used for w.dens
#' @param gamma_b the rate values of the gamma function used for w.dens

plot.simOutbreak <- function(R0=2,gamma_a=4,gamma_b=1){
  temp.sim <- simOutbreak(R0,dgamma(1:10,gamma_a,gamma_b))
  if(length(temp.sim$id)==1) return("No outbreak")
  temp.net <- graph.data.frame(cbind(temp.sim$ances[-1],temp.sim$id[-1]),directed=T)
  return(plot(temp.net))
}

simContact <- function(temp.sim,eps=1,avg.contact=10,N=1000){
  
  generate.contact <- function(id){
    num.contact <- rbinom(1,avg.contact*2,0.5)

    certain.contact <- which(temp.sim$ances==id)
    
    potent.contact <- (1:N)[-c(id,certain.contact)]
    
    if(length(certain.contact)>num.contact) num.contact <- length(certain.contact)
    
    potent.contact <- sample(potent.contact,num.contact-length(certain.contact))
    
    if(num.contact==0) return(NULL)
    
    return(matrix(c(rep(id,num.contact),certain.contact,potent.contact),nrow=num.contact))
  }

  true.contact <- lapply(temp.sim$id,generate.contact)
  
  true.contact <- as.data.frame(do.call(rbind,true.contact))
  
  if(nrow(true.contact)==0) return("No contacts")
  
  colnames(true.contact) = c("id","contact")
  
  report.contact <- true.contact[sample(nrow(true.contact),eps*nrow(true.contact)),]
  
  return(report.contact)
}

#Returns the genetic pseudo-likelihood using the number of genetic differences, length of comparable genome,
#number of generations and mutation rate
mu.func <- function(d,l,k,u){
  return(u^(d)*(1-u)^(k*l - d))
}

#Plots the likelihood given fixed d, l and k across a range of values for u, returning most likely u value
plot.mu.func <- function(d,l,k){
temp.seq = seq(0,1,0.01)
temp.output <- mu.func(d,l,k,temp.seq)
plot(temp.seq,temp.output,type='l',xlab="Mutation rate",ylab="Likelihood")
return(temp.seq[which.max(temp.output)])
}

#Returns the probability of infection from simOutbreak given R0 and n (summing over w(t-ti))
p.inf <- function(R0,n){
  return(1 - exp(1)^(-R0/n))
}

#' @param cij is a binary (0,1) describing if contact tracing is observed between i and j
#' @param eps is the contact tracing coverage

omega.3 <- function(cij,eps){
  return(eps^cij + eps*(cij-1))
}