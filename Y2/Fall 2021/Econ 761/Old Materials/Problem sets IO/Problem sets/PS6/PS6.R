#################
# Problem set 6 #
# Natalia Serna #
#################

rm(list=ls())
library(data.table)
library(MASS)
library(doParallel)
library(oglmx)
library(Matrix)
library(margins)
library(ggplot2)
cl <- makeCluster(2)
registerDoParallel(cl)
setwd("~/Documents/Ph.D/Econ 761 - Industrial Organization Theory/Problem sets/PS6/")
air <- read.csv("airlines.csv", sep=",", header=F)
colnames(air) <- c("citypair", "yr_q", "car", "pax", "pot96", "pot97", "enter", "city2",
                   "city972", "distance", "tourist", "basepop", "refpop", "paxtot", "pop",
                   "enter96", "totenter", "cityN2", "cityN1", "numroute", "sharepax", "sharepaxdist",
                   "totpotential", "herfCityPair", "incumbents", "totsinglepot", "dist2")
air <- data.table(air)

##### Part 1. Summary statistics
air[, mean(enter), by=yr_q]
air[, mean(enter), by=list(yr_q, city2)]
air[, mean((totenter)/cityN2), by=yr_q]
air[, list(sum(numroute), sum(pax)), by=list(yr_q, enter, enter96)]
air[, list(mean(pop), mean(distance), mean(tourist), mean(city2), mean(sharepaxdist), mean(numroute))]
air[, list(sd(pop), sd(distance), sd(tourist), sd(city2), sd(sharepaxdist), sd(numroute))]
air[, list(mean(pop), mean(herfCityPair), mean(tourist), mean(incumbents)), by=list(enter, enter96,yr_q)]
a <- air[which(enter==1),incumbents]
b <- air[which(enter==0),incumbents]
t.test(a,b)

##### Part 2.
# Select subsample used by Berry (1992): all potential entrants with presence on at least one the end points
# and all incumbents
y <- air[which(yr_q==19972 & pot96==1), "enter"]
x1 <- air[which(yr_q==19962 & pot96==1), c("city2", "sharepaxdist", "numroute", "sharepax")]
x2 <- air[which(yr_q==19972 & pot96==1), c("pop", "distance", "tourist")]
dat <- data.frame(cbind(y,x1,x2))
dat$distance2 <- (dat$distance)^2
probit1 <- glm(enter~pop+distance+distance2+tourist+city2+sharepaxdist, 
               family = binomial(link = "probit"), data=dat)
probit2 <- glm(enter~pop+distance+distance2+tourist+city2+numroute, 
               family = binomial(link = "probit"), data=dat)
probit3 <- glm(enter~pop+distance+city2+sharepaxdist, 
               family = binomial(link = "probit"), data=dat)
-2*logLik(probit1)
-2*logLik(probit2)
-2*logLik(probit3)
mfx <- marginal_effects(probit3)
colMeans(mfx)

##### Part 3.
# Column 1 of table 6 in Berry (1992) has each observation being a market. So we need to
# collapse the information by markets
x <- air[, list(mean(pop), mean(distance), mean(totenter), mean(totpotential)), by=c("yr_q", "citypair")]
setnames(x, c("yr_q", "citypair", "pop", "distance", "totenter", "totpotential"))
y <- x[which(yr_q==19972), "totpotential"]
x <- x[which(yr_q==19972), c("pop", "distance", "totenter")]
x[, ltotenter:=log(totenter)]
dat <- cbind(y,x)
op <- oglmx(totpotential~pop+distance+ltotenter, data=dat, link="probit", 
                 constantMEAN=T, constantSD=F, delta=0, threshparam=NULL)
summary(op)

# Testing the restriction of additional entrants
mfxo <- margins.oglmx(op)
mf <- NULL
for(i in 1:length(mfxo)){
  m <- mfxo[[i]]
  mf <- rbind(mf, cbind(m[3,1], m[3,2]))
}
mf <- data.frame(seq(9,21,1), mf)
colnames(mf) <- c("n", "MF", "se")
ggplot(mf, aes(x=n, y=MF)) + geom_errorbar(aes(ymin=MF-1.96*se, ymax=MF+1.96*se), width=.1) +
  geom_line() + geom_point() + theme_bw()

##### Part 4. 
# Second column of table 6
source("oprobit.R")
source("vcovoprobit.R")
air[, constant:=1]
x <- air[, c("yr_q", "citypair", "constant", "pop", "distance", "city2", "sharepaxdist", "totenter", "totpotential")]
y <- x[which(yr_q==19972), c("citypair", "totpotential")]
y <- y[, mean(totpotential), by="citypair"]
setnames(y, c("citypair", "totpotential"))
x1 <- x[which(yr_q==19972), c("citypair", "constant", "pop", "distance","totenter")]
x1 <- x1[, list(mean(constant), mean(pop), mean(distance), mean(totenter)), by="citypair"]
setnames(x1, c("citypair", "constant", "pop", "distance", "totenter"))
x2 <- x[which(yr_q==19962), c("citypair","city2", "sharepaxdist")]
x1[, ltotenter:=log(totenter)]
x1[, totenter:=NULL]
x <- merge(x2, x1, by="citypair")
x <- merge(x,y, by="citypair")
x <- x[, c("citypair", "constant", "pop", "distance", "city2", "sharepaxdist", "ltotenter", "totpotential")]
x <- data.frame(x)
theta0 <- c(-1,1.2,1.5,3.2,3.2,-0.5)
ophet <- optim(theta0, fn=oprobit, method="Nelder-Mead", x=x, control=list(trace=T))
theta <- ophet$par
save(theta, file="thetaoprobit.Rdata")
vcovphet <- vcovoprobit(theta=theta, x=x)
sqrt(diag(vcovphet))

# Testing the restriction of additional entrants
test <- air[which(yr_q==19972), c("citypair", "city2", "sharepaxdist", "totenter", "enter")]
test[, za:=city2*theta[4]+sharepaxdist*theta[5]]
test[, dn:=theta[6]*(log(totenter)-log(totenter+1))]
xt <- split(test, as.factor(test[,citypair]))
acc <- foreach(i=1:length(xt), .combine=rbind, .packages=c("doParallel"))%dopar%{
  xtt <- xt[[i]]
  xtt <- data.frame(xtt)
  dl <- mean(xtt[,"dn"])
  enter <- xtt[which(xtt$enter==1), c("za")]
  notenter <- xtt[which(xtt$enter==0), c("za")]
  r <- foreach(j=1:length(notenter), .combine='+')%dopar%{
    ri <- sum(as.numeric(enter-notenter[j]>dl))
    ri
  }
  c(sum(r), length(enter)*length(notenter))
}
rej <- 1-(sum(acc[,1])/sum(acc[,2]))

##### Part 5.
source("entrysimul.R")
source("combinesim.R")
air[, constant:=1]
x <- air[, c("yr_q", "citypair", "constant", "pop", "distance", "city2", "sharepaxdist", 
             "totenter", "enter", "herfCityPair", "totpotential", "totsinglepot")]
y <- x[which(yr_q==19972), c("citypair", "totenter")]
y <- y[, mean(totenter), by="citypair"]
setnames(y, c("citypair", "totenter"))
enter <- x[which(yr_q==19972), c("enter")]
x <- x[which(yr_q==19962), c("citypair", "constant", "pop", "distance", "city2", "sharepaxdist", 
                             "totenter", "herfCityPair", "totpotential", "totsinglepot")]
x[, ltotenter:=log(totenter)]
x[, totenter:=NULL]
x <- cbind(x, enter)
x <- merge(x,y, by="citypair")
pot <- x[, mean(totpotential), by="citypair"]
pot <- data.frame(pot[,"V1"])
pot <- pot[,1]
z <- data.table(cbind(x[, "herfCityPair"]/10000, log(x[, "totpotential"]), log(x[,"totsinglepot"])))
setnames(z, c("herfindahl", "totpotential", "totsinglepot"))
x[, herfCityPair:=NULL]
x[, totpotential:=NULL]
x[, totsinglepot:=NULL]
x <- data.frame(x)
theta <- c(0.2,0.2,0.2,0.2,0.2,-0.5,0.3)

entrys <- optim(theta, fn=entrysimul, method="Nelder-Mead", x=x, z=z, pot=pot, R=20, order=T, A=diag(30),
                combinesim=combinesim, control=list(maxit=1))

gb1 <- foreach(i=1:ncol(pa), .combine=cbind, .packages=c("doParallel"))%dopar%{
  gbi <- pa[,i]
  gbi <- sweep(x1a, 1, gbi, `*`)
  gbi
}
gb2 <- sweep(x2a, 1, v0a, `*`)

gba <- cbind(rep(gb2[,1], each=31), rep(gb2[,2], each=31), rep(gb2[,3], each=31),
             rep(gb2[,4], each=31), rep(gb2[,5], each=31), rep(gb2[,6], each=31),
             gb1)

A <- (t(gba)%*%gba)/184
A <- solve(A)

theta <- c(-1,1.2,1.5,3.2,3.2,-0.5, 0.1)
entryf <- optim(theta, fn=entrysimul, method="Nelder-Mead", x=x, z=z, pot=pot, R=20, order=T, A=A,
                combinesim=combinesim, control=list(trace=T)) 
theta <- entryf$par
save(theta, file="thetasimul.Rdata")
stopCluster(cl)
