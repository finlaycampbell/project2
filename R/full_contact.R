

#################################################################################################
######### Finlay's ongoing contribution to Outbreaker2, working on incoporating contact ######### 
######### tracing data into the inference of posterior transmission trees.              ######### 
#################################################################################################

##################################################
######### LOADING LIBRARIES AND DATASETS ######### 
##################################################

library(visNetwork)
library(ggplot2)
library(reshape2)
devtools::load_all()

######################################
######### DEFINING FUNCTIONS #########
######################################

## A function to analyse the accuracy and precision of outbreaker inference 
result.analysis <- function(result,true.outbreak,plot=TRUE,print=FALSE){
  id <- seq_len(true.outbreak$n)
  adder <- which(names(result)=="alpha.1")-1
  samples <- length(result$step)
  
  #Determine the modal transmission network
  network <- data.frame(from=do.call(rbind,lapply(id, function(i) ({
    modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
    if(length(modal.ances)==0) return(NA) else return(modal.ances)
  }))), to=id) 
  
  
  #Plot proposed transmission network
  
  import <- which(is.na(network$from))
  
  #Define the indices to call the times of infection
  t.inf.index <- which(names(result)=="t.inf.1")-1+id
  
  #Determine the median posterior time of onset for plotting
  onset <- unlist(lapply(result[t.inf.index],median))
  
  #Scale onset to begin at 0
  onset <- onset - min(onset)
  
  nodes <- data.frame(id=id,label=id)
  nodes$color <- gray(1/3+2*onset/(3*max(onset)))
  nodes$color[import] <- "#e60000"
  
  plot.network <- visNetwork(nodes,network) %>%
    visNodes(shape="ellipse",font=list("size"=25),borderWidth=2) %>%
    visEdges(arrows="to",color="#e60000")
    
  if(plot) print(plot.network)
  
  #Determine confidence in our results

  transmission.id <- id[-sapply(result[id+adder],function(i) any(is.na(i)))]

  confidence <- round(mean(sapply(transmission.id,function(i) mean(result[[i+adder]]==true.outbreak$ances[i]))),2)
  
  entropy <- round(mean(sapply(result[transmission.id+adder],function(i) {fk <- table(i)/sum(table(i))
                                                                                 -sum(log(fk)*fk)})),2)
   
  #Determine the proportion of correctly inferred ancestries
  num.correct <-  sum(true.outbreak$ances==network$from,na.rm=TRUE)
  num.correct <- num.correct + sum(is.na(true.outbreak$ances[is.na(network$from)]))
  accuracy <- round(num.correct/nrow(network),2)
  
  if(print) print(paste(100*accuracy,"% of transmission tree correctly inferred | Average confidence: ",confidence,sep=""))
  
  out <- list(analysis=data.frame(accuracy=accuracy,confidence=confidence,entropy=entropy),plot.network=plot.network)
  
  return(out)
}

## A function to simulate contact tracing data (CTD)
simCTD <- function(temp.outbreak,eps=1,xi=0,plot=FALSE,print.ratio=FALSE){
  
  if(temp.outbreak$n==1) stop("No transmission observed")
  
  is.contact <- function(pair) return(any(pair[1] == infec.contact[,1] & pair[2] == infec.contact[,2]))
  
  accept.reject <- function(pair){
    if(pair[3]) return(runif(1,0,1) < eps)
    else return(runif(1,0,1) < xi*eps)
  }
  
  import <- which(is.na(temp.outbreak$ances))
  infec.contact <- cbind(temp.outbreak$ances[-import],temp.outbreak$id[-import])
  
  potent.CTD <- as.data.frame(t(combn(temp.outbreak$id,2)))
  colnames(potent.CTD) = c("i","j")
  potent.CTD <- cbind(potent.CTD,contact=apply(potent.CTD,1,is.contact))
  potent.CTD <- cbind(potent.CTD,accept=apply(potent.CTD,1,accept.reject))
  
  CTD <- potent.CTD[potent.CTD$accept,1:2]
  rownames(CTD) <- NULL

  if(plot){
    plot.CTD <- potent.CTD[potent.CTD$accept | potent.CTD$contact,]
    
    nodes.color <- gray(1/3+2*temp.outbreak$onset/(3*tail(temp.outbreak$onset,1)))
    nodes.color[import] <- "#e60000"
    
    nodes.label = temp.outbreak$id
    
    edges.color <- rep("#e60000",nrow(plot.CTD))
    edges.color[!plot.CTD$contact] <- "blue"
    edges.arrows <- rep("to",nrow(plot.CTD))
    edges.arrows[!plot.CTD$contact] <- FALSE
    edges.dashes <- plot.CTD$contact & !plot.CTD$accept
    edges.width <- rep(3,nrow(plot.CTD))
    edges.width[!plot.CTD$accept]  <- 0.5
    
    nodes <- data.frame(id=temp.outbreak$id,label=nodes.label,color=nodes.color)
    edges <- data.frame(from=plot.CTD$i,to=plot.CTD$j,dashes=edges.dashes,arrows=edges.arrows,
                        color=edges.color,width=edges.width)
    network <- visNetwork(nodes,edges,main="Simulated Contact Network") %>%
      visNodes(shape="ellipse",font=list("size"=25),borderWidth=2)
    print(network)
  }
  
  return(CTD)
}

## A function to compare accuracy and precision of outbreaker and CTD.outbreaker
compare.outbreakers <- function(runs=50,min.cases=10,n.hosts=15,eps=0.8,xi=0.01,w.dens=dgamma(1:20,2,0.05),
                                 R0=2,mu.transi=1e-05,rate.import.case=0,print=FALSE,seq.length=10000,no.dna=FALSE){

  analysis <- data.frame(matrix(nrow=runs,ncol=8))
  colnames(analysis) <- c("accuracy","CTD.accuracy","confidence","CTD.confidence",
                          "entropy","CTD.entropy","time","CTD.time")
  
  store.outbreak <- store.CTD <- store.their.result <- store.our.result <- list()
  
  counter <- 1

  while(counter<=runs){
    
    outbreak <- outbreaker::simOutbreak(R0,w.dens,n.hosts=n.hosts,mu.transi=mu.transi,
                                        rate.import.case=rate.import.case,seq.length=seq.length)
    if(outbreak$n<min.cases) next
    
    CTD <- simCTD(outbreak,plot=FALSE,eps=eps,xi=xi)
    
    their.dna <- outbreak$dna
    
    if(no.dna) their.dna <- NULL 
    
    their.time <- system.time(
      their.result <- outbreaker(data=list(dates=outbreak$onset,w.dens=w.dens,dna=their.dna),
                                 config=list(n.iter=2e4,sample.every=200)))
    
    their.analysis <- result.analysis(their.result,outbreak,print=FALSE,plot=FALSE)
    analysis[counter,c(1,3,5,7)] <- c(their.analysis$analysis,their.time[1])
    
    our.time <- system.time(
      our.result <- outbreaker(data=list(dates=outbreak$onset,dna=outbreak$dna,w.dens=w.dens,CTD=CTD),
                               config=list(n.iter=2e4, sample.every=200)))
    
    our.analysis <- result.analysis(our.result,outbreak,print=FALSE,plot=FALSE)
    analysis[counter,c(2,4,6,8)] <- c(as.numeric(our.analysis$analysis),our.time[1])
    
    store.outbreak[[counter]] <- outbreak
    store.CTD[[counter]] <- CTD    
    store.their.result[[counter]] <- their.result
    store.our.result[[counter]] <- our.result
    
    counter <- counter + 1
    
  }
  
  plot.out <- analysis
  plot.out[,7:8] <- analysis[,7:8]/max(analysis[,7:8])
  plot.out <- as.data.frame(reshape2::melt(plot.out))
  
  plot.out$variable <- c(rep("Accuracy",2*runs),rep("Confidence",2*runs),rep("Entropy",2*runs),rep("Time",2*runs))
  plot.out$Model <- factor(rep(c(rep("Original",runs),rep("CTD",runs)),4),levels=c("Original","CTD"))
  
  p <- ggplot(plot.out) + geom_violin(aes(x=variable,y=value,fill=Model),scale="width") + ylim(0.,1) +
    ylab("Value") + theme_set(theme_gray(base_size = 18)) + theme(axis.title.x=element_blank())
  
  param <- list(n.hosts=n.hosts,min.cases=min.cases,R0=R0,w.dens=w.dens,n.hosts=n.hosts,
                mu.transi=mu.transi,rate.import.case=rate.import.case,eps=eps,xi=xi,runs=runs)
  
  out <- list()
  out$param <- param
  out$outbreak <- store.outbreak
  out$CTD <- store.CTD
  out$their.result <- store.their.result
  out$our.result <- store.our.result
  out$analysis <- analysis
  out$plot <- p
  
  return(out)
}

## A function to compare accuracy and precision of CTD.outbreaker with DNA and with no DNA
compare.outbreakers2 <- function(runs=50,min.cases=10,n.hosts=15,eps=0.8,xi=0.01,w.dens=dgamma(1:20,2,0.05),
                                R0=2,mu.transi=1e-05,rate.import.case=0,print=FALSE,seq.length=10000,no.dna=FALSE){
  
  analysis <- data.frame(matrix(nrow=runs,ncol=8))
  colnames(analysis) <- c("tc.accuracy","tcg.accuracy","tc.confidence","tcg.confidence",
                          "tc.entropy","tcg.entropy","tc.time","tcg.time")
  
  store.outbreak <- store.CTD <- store.tc.result <- store.tcg.result <- list()
  
  counter <- 1
  
  while(counter<=runs){
    
    outbreak <- outbreaker::simOutbreak(R0,w.dens,n.hosts=n.hosts,mu.transi=mu.transi,
                                        rate.import.case=rate.import.case,seq.length=seq.length)
    if(outbreak$n<min.cases) next
    
    CTD <- simCTD(outbreak,plot=FALSE,eps=eps,xi=xi)

    tc.time <- system.time(
      tc.result <- outbreaker(data=list(dates=outbreak$onset,w.dens=w.dens,dna=NULL,CTD=CTD),
                              config=list(n.iter=2e4,sample.every=200)))
    
    tc.analysis <- result.analysis(tc.result,outbreak,print=FALSE,plot=FALSE)
    analysis[counter,c(1,3,5,7)] <- c(tc.analysis$analysis,tc.time[1])
    
    tcg.time <- system.time(
      tcg.result <- outbreaker(data=list(dates=outbreak$onset,w.dens=w.dens,dna=outbreak$dna,CTD=CTD),
                               config=list(n.iter=2e4, sample.every=200)))
    
    tcg.analysis <- result.analysis(tcg.result,outbreak,print=FALSE,plot=FALSE)
    analysis[counter,c(2,4,6,8)] <- c(as.numeric(tcg.analysis$analysis),tcg.time[1])
    
    store.outbreak[[counter]] <- outbreak
    store.CTD[[counter]] <- CTD    
    store.tc.result[[counter]] <- tc.result
    store.tcg.result[[counter]] <- tcg.result
    
    counter <- counter + 1
    
  }
  
  plot.out <- analysis
  plot.out[,7:8] <- analysis[,7:8]/max(analysis[,7:8])
  plot.out <- as.data.frame(reshape2::melt(plot.out))
  
  plot.out$variable <- c(rep("Accuracy",2*runs),rep("Confidence",2*runs),rep("Entropy",2*runs),rep("Time",2*runs))
  plot.out$Model <- factor(rep(c(rep("Original",runs),rep("CTD",runs)),4),levels=c("Original","CTD"))
  
  p <- ggplot(plot.out) + geom_violin(aes(x=variable,y=value,fill=Model),scale="width") + ylim(0.,1) +
    ylab("Value") + theme_set(theme_gray(base_size = 18)) + theme(axis.title.x=element_blank())
  
  param <- list(n.hosts=n.hosts,min.cases=min.cases,R0=R0,w.dens=w.dens,n.hosts=n.hosts,
                mu.transi=mu.transi,rate.import.case=rate.import.case,eps=eps,xi=xi,runs=runs)
  
  out <- list()
  out$param <- param
  out$outbreak <- store.outbreak
  out$CTD <- store.CTD
  out$tc.result <- store.tc.result
  out$tcg.result <- store.tcg.result
  out$analysis <- analysis
  out$plot <- p
  
  return(out)
}

## A function to compare the relative perfomanrce of outbreaker and CTD.outbreaker over parameter ranges
compare.lineplot <- function(storage,var,x.lab="",...){
  
  iterations <- length(storage)
  runs <- storage[[names(storage)[1]]]$param$runs
  
  relative.analysis <- matrix(ncol=5,nrow=iterations)
  absolute.analysis <- matrix(ncol=9,nrow=iterations)
  full.abs.analysis <- matrix(ncol=9,nrow=iterations*runs)
  
  counter <- 1
  
  for(run in names(storage)){
    analysis <- storage[[run]]$analysis
    mean.analysis <- sapply(analysis,mean)
    
    absolute.analysis[counter,1:8] <- mean.analysis
    full.abs.analysis[(1:runs)+runs*(counter-1),1:8] <- as.matrix(analysis)
    
    for(i in 1:4) relative.analysis[counter,i] <- (mean.analysis[2+2*(i-1)]/mean.analysis[1+2*(i-1)]) - 1
    relative.analysis[counter,5] <- absolute.analysis[counter,9] <- storage[[run]]$param[[var]]
    full.abs.analysis[(1:runs)+runs*(counter-1),9] <- storage[[run]]$param[[var]]
    
    counter <- counter + 1
  }
  
  colnames(relative.analysis) <- c("Accuracy","Entropy","Confidence","Time",var)
  ggplot.rel.analysis <- melt(as.data.frame(relative.analysis),id=var)
  names(ggplot.rel.analysis) <- c(var,"Variable","Value")
  ggplot.rel.analysis$Variable <- factor(ggplot.rel.analysis$Variable,levels=c("Accuracy","Confidence","Entropy","Time"))
  ggplot.rel.analysis$Value <- ggplot.rel.analysis$Value*100
  
  colnames(absolute.analysis) <- c(names(mean.analysis),var)
  ggplot.abs.analysis <- melt(as.data.frame(absolute.analysis),id=var)
  ggplot.abs.analysis$variable <- c(rep("Accuracy",iterations*2),rep("Entropy",iterations*2),
                                 rep("Confidence",iterations*2),rep("Time",iterations*2))
  ggplot.abs.analysis$Model <- factor(rep(c(rep("Original",iterations),rep("CTD",iterations)),4),levels=c("Original","CTD"))
  names(ggplot.abs.analysis) <- c(var,"Variable","Value","Model")
  
  colnames(full.abs.analysis) <- c(names(mean.analysis),var)
  ggplot.full.abs.analysis <- melt(as.data.frame(full.abs.analysis),id=var)
  ggplot.full.abs.analysis$variable <- factor(c(rep("Accuracy",runs*iterations*2),rep("Entropy",runs*iterations*2),
                                    rep("Confidence",runs*iterations*2),rep("Time",runs*iterations*2)),
                                    levels=c("Accuracy","Confidence","Entropy","Time"))
  ggplot.full.abs.analysis$Model <- factor(rep(c(rep("Original",runs*iterations),rep("CTD",runs*iterations)),4),levels=c("Original","CTD"))
  names(ggplot.full.abs.analysis) <- c(var,"Variable","Value","Model")
  
  
  p <- ggplot(ggplot.rel.analysis,aes(x=ggplot.rel.analysis[[var]],y=Value,col=Variable)) + geom_line(size=1.5) + 
       xlab(x.lab) + ylab("Relative change to Outbreaker (%)")
  p
  
  q <- ggplot(ggplot.abs.analysis,aes(x=ggplot.abs.analysis[[var]],y=Value,col=Model)) + geom_line(size=1.5) + 
       facet_grid(Variable ~ .,scales="free") + xlab(x.lab)
  q
  
  r <- ggplot(ggplot.full.abs.analysis,aes(x=ggplot.full.abs.analysis[[var]],y=Value,col=Model)) + stat_smooth() +
       facet_grid(Variable ~ .,scales="free") + xlab(x.lab)
  r
  
  #ggsave("C:\\Users\\Finlay\\Documents\\Projects\\project2\\figs\\relative.analysis.png",p)
  #ggsave("C:\\Users\\Finlay\\Documents\\Projects\\project2\\figs\\absolute.analysis.png",q)
  #ggsave("C:\\Users\\Finlay\\Documents\\Projects\\project2\\figs\\smooth.analysis.png",r)
  
  return(list("rel.plot"=p,"abs.plot"=q,"smooth.plot"=r))
}

## This function will run compare.outbreakers over the parameter range provided in ... in the form var=c(0.1,0.2,...)
## It returns a list storing all parameters, outbreaks, outbreaker results, analyses and plots 
general.comparison <- function(runs,) {

  #Call the list of arguments provided to run.comparison
  args <- as.list(match.call())
  
  #Remove the first element from this list (not used as argument input for next function)
  args[[1]] <- NULL
  
  #Define our variable as the first element in ...
  var <- names(list(...))[1]
  
  #Extract the values var takes
  seq <- list(...)[[var]]

  for(value in seq){
    arg.mod <- args
    arg.mod[[var]] <- value
    
    print(paste(round(which(seq==value)/length(seq)*100),"%",sep=""))

    storage[[paste(var,round(value*100,1),sep="")]] <- do.call(compare.outbreakers,arg.mod) 
  }
    
  plots <- compare.lineplot(storage,var)
  
  out <- list(data=storage,plots=plots)

}

grid.comparison <- function() {
  
  runs <- 50
  n.hosts <- 200
  min.cases <- 10
  rate.import.case <- 0
  
  ebola.w.dens <- 
  ebola.R0 <- 1.5
  ebola.mu.transi <- 
  ebola.seq.length <- 
    
  sars.w.dens <- 
  sars.R0 <- 3.5
  sars.mu.transi <- 
  sars.seq.length <- 
    
  eps.seq <- seq(0,1,0.1)  
  xi.seq <- seq(0,1,0.1)
  
  ebola <- list()
  sars <- list()

  for(eps in eps.seq){
    for(xi in xi.seq){
      
      #Avoid running eps=0 and various combinations of xi as they will be identical - run only eps=0/xi=0 once
      if(eps==0 & xi != 0) next
      
      ebola[[paste("eps",eps,"xi",xi,sep="")]] <- compare.outbreakers(runs=runs,min.cases=min.cases,n.hosts=n.hosts,
                                                                      eps=eps,xi=xi,w.dens=ebola.w.dens,R0=ebola.R0,
                                                                      mu.transi=ebola.mu.transi,seq.length=ebola.seq.length)
      
      sars[[paste("eps",eps,"xi",xi,sep="")]] <- compare.outbreakers(runs=runs,min.cases=min.cases,n.hosts=n.hosts,
                                                                      eps=eps,xi=xi,w.dens=sars.w.dens,R0=sars.R0,
                                                                      mu.transi=sars.mu.transi,seq.length=sars.seq.length)
    }
  }
  
  out <- list(ebola=ebola,sars=sars)
  
}


