entryloglik <- function(theta, x, y, rho, R, n){
  
  source("entrymle.R")
  xt <- split(x, as.factor(x$citypair))
  pr <- foreach(t=1:length(xt), .combine=rbind)%dopar%{
    xtt <- xt[[t]]
    xtt <- xtt[,2:ncol(xtt)]
    n <- entrymle(theta=theta, x=xtt, y=yt, rho=rho, R=R, n=nrow(xtt))
    yt <- y[t]
    if(yt==0) yt <- c(1,0,0)
    if(yt==1) yt <- c(0,1,0)
    if(yt==2) yt <- c(0,0,1)
    n <- n^yt
    n <- prod(n)
  }
  pr <- prod(pr)
  pr <- log(pr)
  return(-pr)
}