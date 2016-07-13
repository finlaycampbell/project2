#' project2 - Functions
#' 
#' 
#' 
#' @author F. Campbell

simCTD <- function(temp.outbreak,eps=1,chi=1,ksi=0,plot=FALSE,print.ratio=FALSE){
  
  if(temp.outbreak$n==1) return("No transmission observed")
  
  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
  
  accept.reject <- function(pair){
    if(pair[3]) return(runif(1,0,1) < chi*eps)
    else return(runif(1,0,1) < ksi*chi*eps)
  }
  
  import <- which(is.na(temp.outbreak$ances))
  infec.contact <- cbind(temp.outbreak$ances[-import],temp.outbreak$id[-import])
  
  potent.CTD <- as.data.frame(t(combn(temp.outbreak$id,2)))
  colnames(potent.CTD) = c("from","to")
  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))
  
  CTD <- potent.CTD[potent.CTD$accept,1:3]
  
  plot.CTD <- potent.CTD[potent.CTD$accept | potent.CTD$contact,]
  
  if(print.ratio) print(paste("un/in:",round(sum(!CTD$contact)/sum(CTD$contact),2)))
  
  if(plot){
    #nodes.color <- gray.colors(temp.outbreak$n,start=0,end=0.9)

    nodes.color <- gray(1/3+2*temp.outbreak$onset/(3*tail(temp.outbreak$onset,1)))
    nodes.color[import] <- "#e60000"

    nodes.label = temp.outbreak$id
    
    edges.color <- rep("#e60000",nrow(plot.CTD))
    edges.color[!plot.CTD$contact] <- "blue"
    edges.arrows <- rep("to",nrow(plot.CTD))
    edges.dashes <- plot.CTD$contact & !plot.CTD$accept
    edges.width <- rep(3,nrow(plot.CTD))
    edges.width[!plot.CTD$accept]  <- 0.5
    
    nodes <- data.frame(id=temp.outbreak$id,label=nodes.label,color=nodes.color)
    edges <- data.frame(from=plot.CTD$from,to=plot.CTD$to,dashes=edges.dashes,arrows=edges.arrows,color=edges.color,width=edges.width)
    network <- visNetwork(nodes,edges,main="Simulated Contact Network") %>%
      visNodes(shape="ellipse",font=list("size"=25),borderWidth=2)
    print(network)
  }
  
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



#' omega.3 returns the likelihood of observing cij given the tracing coverage eps
#' @param cij is a binary (0,1) describing if contact tracing is observed between i and j
#' @param eps is the contact tracing coverage

omega.3 <- function(cij,eps){
  return(eps^cij + eps*(cij-1))
}

omega.3.cij <- function(actual.CTD,prop.cij,eps){
    cij <- which(actual.CTD$from==prop.cij[1] & actual.CTD$to==prop.cij[2])
    cij <- length(cij)
    
    if(cij > 1) {return(print("Double contact reported"));return(NA)}
    
    return(eps^cij + eps*(cij-1))
  }

#' ll.CTD returns the pseudo-likelihood of proposed CTD given the actual, observed CTD
ll.CTD <- function(actual.CTD,prop.contact,eps){
  
  temp.ll = 1
  
  for(i in 1:nrow(prop.contact)){
    temp.ll <- temp.ll*omega.3.cij(actual.CTD,prop.contact[i,],eps)
  }
  
  return(temp.ll)
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