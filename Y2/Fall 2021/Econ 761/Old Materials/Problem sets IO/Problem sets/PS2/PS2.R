###################
## Problem set 2 ##
## Natalia Serna ##
###################

rm(list=ls())
setwd("~/Documents/Ph.D/Econ 761 - Industrial Organization Theory/Problem sets/PS2/")
set.seed(111)

#Part 2.a and 2.b

#Data generating process
c0=1
c1=0.9
xi=0
f=1
b0=1
b1=0
eta=0
cities=1000
N=ceiling(10*runif(cities))
collusion=as.numeric(c((N[1:500]<=8), rep(FALSE, 500)))

LernerCN=c1/N
HHI=1/N
Elasticity=rep(1/c1,cities)
LernerM=c1
Lerner=collusion*LernerM+(1-collusion)*LernerCN
oblerner=log(Lerner)+0.1*(runif(cities)-0.5)
loghhi=log(HHI)

#Part 2.c
data = data.frame(oblerner, loghhi, collusion)

#Regressions
col = lm(oblerner~loghhi, data=data[1:500,])
cou = lm(oblerner~loghhi, data=data[501:1000,])
all = lm(oblerner~loghhi, data=data)

#Part 2.c
#Hypothesis testing
#Collusion
sumcol=summary(col)$coefficients
tcol=abs(sumcol[2,1]-1)/sumcol[2,2]
rejectcol=(tcol>=qnorm(0.95))
#Cournot
sumcou=summary(cou)$coefficients
tcou=abs(sumcou[2,1]-1)/sumcou[2,2]
rejectcou=(tcou>=qnorm(0.95))
#All
sumall=summary(all)$coefficients
tall=abs(sumall[2,1]-1)/sumall[2,2]
rejectall=(tall>=qnorm(0.95))

#Part 2.d
a0=3
a1=1
v=0
source("functions.R")

Lerner=NULL
hhi=NULL
Elasticity=NULL

for(i in 1:cities){
  if(collusion[i]==T){
    Q=(a0+v-b0-eta)/(2*a1)
    q=Q/N[i]
    q=rep(1, N[i])
    p=0.5*(a0+v+b0+eta)
    l1=lerner(p=p, mc=b0+eta)
    Lerner[i]=l1
    h=HHI(q=q, Q=Q)
    hhi[i]=h
    Elasticity=rbind(Elasticity,(1/a1)*((a0-a1*Q+v)/Q))
  }else{
    q=(a0+v-b0-eta)/(a1*(N[i]+1))
    Q=N[i]*q
    q=rep(q, N[i])
    p=(a0+v+N[i]*(b0+eta))/(N[i]+1)
    l1=lerner(p=p, mc=b0+eta)
    Lerner[i]=l1
    h=HHI(q=q, Q=Q)
    hhi[i]=h
    Elasticity=rbind(Elasticity,abs((-1/a1)*((a0-a1*Q+v)/Q)))
  }
}

oblerner=log(Lerner)+runif(cities,-0.05, 0.05)
loghhi=log(hhi)
data = data.frame(oblerner, loghhi, collusion)

#Regressions
col = lm(oblerner~loghhi, data=data[1:500,])
cou = lm(oblerner~loghhi, data=data[501:1000,])
all = lm(oblerner~loghhi, data=data)

#Part 2.c
#Hypothesis testing
#Collusion
sumcol=summary(col)$coefficients
tcol=abs(sumcol[2,1]-1)/sumcol[2,2]
rejectcol=(tcol>=qnorm(0.95))
#Cournot
sumcou=summary(cou)$coefficients
tcou=abs(sumcou[2,1]-1)/sumcou[2,2]
rejectcou=(tcou>=qnorm(0.95))
#All
sumall=summary(all)$coefficients
tall=abs(sumall[2,1]-1)/sumall[2,2]
rejectall=(tall>=qnorm(0.95))

#Part 3.a

a0=5
a1=1
f=1
b0=1
b1=0
cities=1000

eta=0
v=runif(cities,-1,1)
N=NULL
for(i in 1:cities){N=rbind(N,((a0+v[i]-b0-eta)/sqrt(a1*f))-1)}

Lerner=NULL
hhi=NULL
Elasticity=NULL
for(i in 1:cities){
  q=(a0+v[i]-b0-eta)/(a1*(N[i]+1))
  Q=N[i]*q
  q=rep(q, N[i])
  p=(a0+v[i]+N[i]*(b0+eta))/(N[i]+1)
  l1=lerner(p=p, mc=b0+eta)
  Lerner=rbind(Lerner,l1)
  h=HHI(q=q, Q=Q)
  hhi=rbind(hhi,h)
  Elasticity=rbind(Elasticity,abs((-1/a1)*((a0-a1*Q+v[i])/Q)) )
}

oblerner=log(Lerner)+runif(cities, -0.05, 0.05)
loghhi=log(hhi)
data = data.frame(oblerner, loghhi)

#Regressions
end_1 = lm(oblerner~loghhi, data=data)

#Part 3.b
eta=runif(cities,-1,1)
v=0
N=NULL
for(i in 1:cities){N=rbind(N,((a0+v-b0-eta[i])/sqrt(a1*f))-1)}

Lerner=NULL
hhi=NULL
Elasticity=NULL
for(i in 1:cities){
  q=(a0+v-b0-eta[i])/(a1*(N[i]+1))
  Q=N[i]*q
  q=rep(q, N[i])
  p=(a0+v+N[i]*(b0+eta[i]))/(N[i]+1)
  l1=lerner(p=p, mc=b0+eta[i])
  Lerner=rbind(Lerner,l1)
  h=HHI(q=q, Q=Q)
  hhi=rbind(hhi,h)
  Elasticity=rbind(Elasticity,abs((-1/a1)*((a0-a1*Q+v)/Q))  )
}

oblerner=log(Lerner)+runif(cities, -0.05, 0.05)
loghhi=log(hhi)
data = data.frame(oblerner, loghhi)

#Regressions
end_2 = lm(oblerner~loghhi, data=data)
