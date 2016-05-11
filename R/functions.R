#' project2 - Functions
#' 
#' 
#' 
#' @author F. Campbell



#' simCTD returns a simulated contact tracing data set, given the known transmission network
#' @param temp.outbreak the known or simulated dataset of type simOutbreak
#' @param eps the contact tracing coverage
#' @param avg.contact the average number of total (infectious and non-infectious) contacts
#' @param N the total number of individuals in the closed system

simCTD <- function(temp.outbreak,eps=1,avg.contact=10,N=1000){
  
  generate.contact <- function(id){
    num.contact <- rbinom(1,avg.contact*2,0.5)

    certain.contact <- which(temp.outbreak$ances==id)
    
    potent.contact <- (1:N)[-c(id,certain.contact)]
    
    if(length(certain.contact)>num.contact) num.contact <- length(certain.contact)
    
    potent.contact <- sample(potent.contact,num.contact-length(certain.contact))
    
    if(num.contact==0) return(NULL)
    
    return(matrix(c(rep(id,num.contact),certain.contact,potent.contact),nrow=num.contact))
  }

  true.contact <- lapply(temp.outbreak$id,generate.contact)
  
  true.contact <- as.data.frame(do.call(rbind,true.contact))
  
  if(nrow(true.contact)==0) return("No contacts")
  
  colnames(true.contact) = c("id","contact")
  
  report.contact <- true.contact[sample(nrow(true.contact),eps*nrow(true.contact)),]
  
  return(report.contact)
}



#' infec.contact.reported function returns the proportion of infectious contacts in temp.outbreak reported by reported.contact
#' @param temp.outbreak the known or simulated dataset of type simOutbreak
#' @param the simulated contact tracing dataset returned by reported.contact

infec.contact.reported <- function(temp.outbreak,temp.CTD){
  
  match.contact <- function(i){
  return(temp.outbreak$id[which(temp.outbreak$ances==i)] %in% subset(temp.CTD$contact,temp.CTD$id==i))
}

  if(temp.outbreak$n<=1) return("No true contacts")
  infect.contact.reported <- do.call(sum,lapply(temp.outbreak$id,match.contact))/(temp.outbreak$n-sum(is.na(temp.outbreak$ances)))
  
  return(infect.contact.reported)
}


#' omega.3 returns the likelihood of observing cij given the tracing coverage eps
#' @param cij is a binary (0,1) describing if contact tracing is observed between i and j
#' @param eps is the contact tracing coverage

omega.3 <- function(cij,eps){
  return(eps^cij + eps*(cij-1))
}



#' plot.simOutbreak plots the simulated transmission network from simOutbreak
#' @param R0 the effective reproduction number
#' @param gamma_a the shape value of the gamma function used for w.dens
#' @param gamma_b the rate values of the gamma function used for w.dens

plot.simOutbreak <- function(R0=2,gamma_a=4,gamma_b=1){
  temp.outbreak <- simOutbreak(R0,dgamma(1:10,gamma_a,gamma_b))
  if(length(temp.outbreak$id)==1) return("No outbreak")
  temp.net <- graph.data.frame(cbind(temp.outbreak$ances[-1],temp.outbreak$id[-1]),directed=T)
  return(plot(temp.net))
}



#' mu.func returns the genetic pseudo-likelihood of a proposed transmission network
#' @param d the number of genetic differences
#' @param l the length of comparable genome
#' @param k the number of generations
#' @param u the mutation rate

mu.func <- function(d,l,k,u){
  return(u^(d)*(1-u)^(k*l - d))
}



#' plot.mu.func plots the likelihood of mu.func given fixed d, l and k across a range of values for u, returning most likely u value
#' @param d the number of genetic differences
#' @param l the length of comparable genome
#' @param k the number of generations

plot.mu.func <- function(d,l,k){
temp.seq = seq(0,1,0.01)
temp.output <- mu.func(d,l,k,temp.seq)
plot(temp.seq,temp.output,type='l',xlab="Mutation rate",ylab="Likelihood")
return(temp.seq[which.max(temp.output)])
}