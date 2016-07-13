source("R/functions.R")

library(visNetwork)
library(outbreaker2)

data <- outbreaker.data(dates=fakeOutbreak$collecDates,dna=fakeOutbreak$dat$dna,w.dens=fakeOutbreak$w)

config <- outbreaker.config()

param <- outbreaker.mcmc.init(data,config)
