vcovoprobit <- function(theta, x){
  
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
    drt <- -dnorm(-vy)+dnorm(-vy1)
    c(y, drt)
  }
  
  prr <- rep(pr[,2], each=nrow(xt[[1]]))
  xt <- x
  xt$dr <- prr
  xt <- split(xt, as.factor(xt$citypair)) 
  
  vcovi <- foreach(i=1:length(xt), .combine=rbind, .packages=c("doParallel"))%dopar%{
    xtt <- xt[[i]]
    xtt$citypair <- NULL
    xtt$totpotential <- NULL
    dr <- data.matrix(xtt$dr)
    xtt$dr <- NULL
    xtt <- data.matrix(xtt)
    vi <- t(dr)%*%xtt
    vi
  }
  
  vcovi <- colSums(vcovi)/(nrow(pr))
  vcov <- vcovi%*%t(vcovi)
  return(vcov)
  
}