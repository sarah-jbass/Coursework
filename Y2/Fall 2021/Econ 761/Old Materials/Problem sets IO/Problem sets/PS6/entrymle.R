entrymle <- function(theta, x, y, rho, R, n){
  
  col <- which(names(x)=="logn")
  sig <- matrix(rho, n, n)
  diag(sig) <- rep(1, n)
  
  n <- foreach(i=1:R, .combine=rbind, .packages=c("doParallel"))%dopar%{
    
    u <- mvrnorm(n=n, mu=rep(0,n), Sigma=sig)  
    x[,col] <- 1
    v <- x%*%t(theta)  
    id0 <- as.numeric(ek<=-v)
    if(sum(id0)==n) n0 <- 1 else n0 <- 0
    
    x[,col] <- 2
    n1 <- foreach(k=1:n, .combine=rbind)%dopar%{
      x[k,col] <- 1
      v <- x%*%t(theta)
      id1 <- as.numeric(ek[k]>-v[k])
      id0 <- as.numeric(ek[-k]<=-v[-k])
      id1 <- id1+sum(id0)
      if(id1==n) ni1 <- 1 else ni1 <- 0
    }
    n1 <- sum(n1)
    c(n0, n1)
    
  }
  
  n <- colMeans(n)
  n <- c(n, 1-sum(n))
  return(n)
}