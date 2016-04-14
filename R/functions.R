#DEFINING FUNCTIONS

say_hello <- function(your_name){
  return(paste("Hello ",your_name,", how are you doing?",sep=""))
}

gamma_plot <- function(gamma_a,gamma_b){
  temp_seq = seq(0,10,0.01)
  return(qplot(temp_seq,dgamma(temp_seq,gamma_a,gamma_b),xlab="Value",ylab="Probability density"))
}

#helloooooooooooooooooo