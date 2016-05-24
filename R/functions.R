#' project2 - Functions
#' 
#' 
#' 
#' @author F. Campbell

simCTD <- function(temp.outbreak,eps=1,chi=1,ksi=0,plot=FALSE,print.ratio=FALSE){
  
  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
  
  accept.reject <- function(pair){
    if(pair[3]) return(runif(1,0,1) < chi*eps)
    else return(runif(1,0,1) < ksi*eps)
  }

  na.contact <- which(is.na(temp.outbreak$ances))
  infec.contact <- cbind(temp.outbreak$ances[-na.contact],temp.outbreak$id[-na.contact])
  
  potent.CTD <- as.data.frame(t(combn(temp.outbreak$id,2)))
  colnames(potent.CTD) = c("from","to")
  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))
  
  CTD <- potent.CTD[potent.CTD$accept,1:3]
  
  if(plot){
    nodes = data.frame(id=CTD$from)
    edges = data.frame(from=CTD$from,to=CTD$to,dashes=!CTD$contact)
    visNetwork(nodes,edges,width="100%")
  }
  
  if(print.ratio) print(paste("un/in:",round(sum(!CTD$contact)/sum(CTD$contact),2)))
  
  return(CTD)
}
 


#' simCTD2 returns a simulated contact tracing data (CTD) set, given the known transmission network
#' @param temp.outbreak the known or simulated dataset of type simOutbreak
#' @param eps the contact tracing coverage
#' @param avg.contact the average number of total (infectious and non-infectious) contacts
#' @param N the total number of individuals in the closed system
#' @param plot a boolean determining if the simulated CTD should be plotted as a network
#' @param infec.reported a boolean determining if the proportion of infectious contacts reported should be printed

simCTD2 <- function(temp.outbreak,eps=1,avg.contact=10,N=1000,plot=FALSE,infec.reported=FALSE){
  
  #generate.CTD returns a matrix of CTD for a single id from temp.outbreak
  generate.CTD <- function(id){
    num.contact <- rbinom(1,avg.contact*2,0.5)

    infec.contact <- which(temp.outbreak$ances==id)
    
    #If num.contact < num.infec.contact, set num.contact = num.infect.contact to force reporting of infectious contacts
    if(num.contact<length(infec.contact)) num.contact <- length(infec.contact)
    
    if(num.contact==0) return(NULL)
    
    potential.non.infec.contact <- (1:N)[-c(id,infec.contact)]
    
    non.infec.contact <- sample(potential.non.infec.contact,num.contact-length(infec.contact))
    
    return(matrix(c(rep(id,num.contact),infec.contact,non.infec.contact),nrow=num.contact))
  }

  full.CTD <- lapply(temp.outbreak$id,generate.CTD)
  
  full.CTD <- as.data.frame(do.call(rbind,full.CTD))
  
  if(nrow(full.CTD)==0) return("No contacts")
  
  colnames(full.CTD) = c("id","contact")
  
  obsv.CTD <- full.CTD[sample(nrow(full.CTD),eps*nrow(full.CTD)),]
  
  if(plot) plot(graph.data.frame(obsv.CTD))
  
  if(infec.reported) print(paste("Infectious contacts reported: ",round(infec.contact.reported(temp.outbreak,obsv.CTD)*100,0),"%",sep=""))
  
  return(obsv.CTD)
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



#' ll.CTD returns the pseudo-likelihood of proposed CTD given the actual, observed CTD
ll.CTD <- function(actual.CTD,prop.contact,eps){
  
  omega.3.cij <- function(actual.CTD,prop.cij,eps){
    cij <- which(actual.CTD$id==prop.cij[1] & actual.CTD$ances==prop.cij[2])
    cij <- length(cij)
    
    if(cij > 1) {return(print("Double contact reported"));return(NA)}
    
    return(eps^cij + eps*(cij-1))
  }
  
  temp.ll = 1
  
  for(i in 1:nrow(prop.contact)){
    temp.ll <- temp.ll*omega.3.cij(actual.CTD,prop.contact[i,],eps)
  }
  
  return(temp.ll)
}



#' plot.simOutbreak plots the simulated transmission network from simOutbreak
#' @param R0 the effective reproduction number
#' @param gamma_a the shape value of the gamma function used for w.dens
#' @param gamma_b the rate values of the gamma function used for w.dens

plot.simOutbreak <- function(temp.outbreak){
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