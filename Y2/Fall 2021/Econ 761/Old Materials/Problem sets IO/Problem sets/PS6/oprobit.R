oprobit <- function(theta, x){
  
  xt <- split(x, as.factor(x$citypair))
  
  pr <- foreach(i=1:length(xt), .combine=rbind, .packages=c("doParallel"))%dopar%{
    xtt <- xt[[i]]
    xtt$citypair <- NULL
    xtt$totpotential <- NULL
    y <- exp(mean(xtt[,ncol(xtt)]))
    xtt <- data.matrix(xtt)
    v <- xtt%*%t(t(theta))
    v <- sort(v, decreasing=T)
    vy <- v[y]
    xtt[,ncol(xtt)] <- log(y+1)
    v <- xtt%*%t(t(theta))
    v <- sort(v, decreasing=T)
    vy1 <- v[y]  
    prt <- (1-pnorm(-vy))-(1-pnorm(-vy1))
    c(y, prt)
  }
  
  loglik <- -2*sum(log(pr[,2]))
  return(loglik)
  
}