library(censReg)
setwd("/Users/arielbuehler/Documents/Documents/Cornell/Research/MonteCarloA")
MPNdata <- read.csv("CrossSectional_RawMPN_3.6.17.csv")
MPNdata$log.adj.mpn <- log10(MPNdata$ADJ.MPN)


### fitdistr
library(fitdistrplus)
MPNdata$logleft <- log10(MPNdata$left)
MPNdata$logright <- log10(MPNdata$right)
MPNdata$logMPN <- log10(MPNdata$MPN)
cens_data <- MPNdata[,c("logleft","logright")]
names(cens_data) <- c("left","right")
fit <- fitdistcens(censdata = cens_data,distr = "norm")
hist(MPNdata$logMPN,freq=F)
curve(dnorm(x,mean=fit$estimate[1],sd = fit$estimate[2]),add=T)
