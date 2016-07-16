## A function to compare accuracy and precision of outbreaker and CTD.outbreaker
compare.outbreakers2 <- function(runs=50,min.cases=10,n.hosts=15,eps=0.8,xi=0.01,w.dens=dgamma(1:20,2,0.05),
                                 R0=2,mu.transi=1e-05,rate.import.case=0,print=FALSE){

  analysis <- data.frame(matrix(nrow=runs,ncol=8))
  colnames(analysis) <- c("accuracy","CTD.accuracy","confidence","CTD.confidence",
                          "entropy","CTD.entropy","time","CTD.time")
  
  store.outbreak <- store.CTD <- store.their.result <- store.our.result <- list()

  counter <- 1
  
  while(counter<=runs){
    
    outbreak <- outbreaker::simOutbreak(R0,w.dens,n.hosts=n.hosts,mu.transi=mu.transi,
                                        rate.import.case=rate.import.case)
    if(outbreak$n<min.cases) next

    CTD <- simCTD(outbreak,plot=FALSE,eps=eps,xi=xi)

    their.time <- system.time(
    their.result <- outbreaker(data=list(dates=outbreak$onset,w.dens=w.dens,dna=outbreak$dna),
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


for(eps in c(0,0.5,1)) out[[paste("eps0",eps*10,sep="")]] <- compare.outbreakers2(runs=3,eps=eps)
