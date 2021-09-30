entrysimul <- function(theta, x, z, R, order=T, A, pot, combinesim){
  
  cat(theta, '\n') 
  flush.console() 
  xt <- split(x, as.factor(x$citypair))
  rho <- theta[length(theta)]
  theta <- theta[1:(length(theta)-1)]
  
  Nh <- foreach(i=1:length(xt), .combine=rbind, .packages=c("doParallel", "base"))%dopar%{
    
    xtt <- xt[[i]]
    xtt$citypair <- NULL
    enter <- xtt$enter
    xtt$enter <- NULL
    xtt <- data.matrix(xtt[,1:(ncol(xtt)-1)])
    n <- nrow(xtt)  
    
    np <- foreach(j=1:R, .combine=combinesim, .packages=c("doParallel", "base"))%dopar%{
      
      u <- rnorm(n+1)
      u0 <- u[1]
      u <- u[-1]
      e <- rho*u0+sqrt(1-rho^2)*u
      xtj <- xtt
      
      nr <- foreach(k=1:pot[i], .combine=rbind)%dopar%{
        xtj[, ncol(xtj)] <- log(k)
        pi <- xtj%*%t(t(theta))+e
        nrj <- sum(as.numeric(pi>=0))
        if(is.na(nrj)) nrj <- 0
        if(nrj>=k) nrj <- k else nrj <- 0
        nrj
      }
      
      nr <- max(nr)
                  
      if(order==T){
        xtj <- xtt
        pi <- xtj%*%t(t(theta))+e
        xtj <- cbind(xtj, seq(1,nrow(xtj)))
        xtj <- cbind(xtj, pi)      
        xtj <- xtj[order(xtj[,ncol(xtj)], decreasing=T),]
        xtj <- cbind(xtj, seq(1, nrow(xtj)))
        xtj <- xtj[order(xtj[,(ncol(xtj)-2)], decreasing=F),]
        ph <- t(t(as.numeric(xtj[,ncol(xtj)]<=nr)))
      }
      
      return(list(nr = nr, ph = ph))
    }
    
    nj <- t(t(rep(np$nr/R, n)))
    p <- np$ph
    p <- enter-p/R
    xtt <- cbind(xtt, seq(1,n))
    xtt <- xtt[order(xtt[,"sharepaxdist"], decreasing=T),]
    xtt <- cbind(xtt, seq(1, n))
    xtt <- xtt[order(xtt[,(ncol(xtt)-1)], decreasing=F),]
    id <- which(xtt[, ncol(xtt)] %in% c(1,2,3,4))
    p1 <- p
    p2 <- p
    p3 <- p
    p4 <- p
    p1[-id[1]] <- 0
    p2[-id[2]] <- 0
    p3[-id[3]] <- 0
    p4[-id[4]] <- 0
    npj <- cbind(nj, p1, p2, p3, p4)
    return(npj)
    
  }
  
  x1 <- x
  x1 <- data.table(x1)
  x1[, citypair:=NULL]
  x1[, constant:=NULL]
  x1[, totenter:=NULL]
  x1[, enter:=NULL]
  x1[, herfindahl:=z[,1]]
  x1 <- data.matrix(x1)
  
  p <- Nh[,2:ncol(Nh)] 
  
  gb1 <- foreach(i=1:ncol(p), .combine=cbind, .packages=c("doParallel"))%dopar%{
    gbi <- p[,i]
    gbi <- sweep(x1, 1, gbi, `*`)
    gbi
  }
  gb1 <- colSums(gb1)/184
  
  v0 <- data.table(cbind(x[, "citypair"], x[,"totenter"]-Nh[,1]))
  setnames(v0, c("citypair", "v0"))
  v0 <- v0[, mean(v0), by="citypair"]
  setnames(v0, c("citypair", "v0"))
  v0 <- data.matrix(v0[,"v0"])
  x2 <- x
  x2 <- data.table(x2)
  x2 <- x2[, herfindahl:=z[,"herfindahl"]]
  x2 <- x2[, totpotential:=z[,"totpotential"]]
  x2 <- x2[, totsinglepot:=z[,"totsinglepot"]]
  x2 <- x2[, list(mean(pop), mean(distance), mean(ltotenter), mean(totpotential), mean(totsinglepot), mean(herfindahl)), by="citypair"]
  setnames(x2, c("citypair", "pop", "distance","ltotenter", "totpotential", "totsinglepot", "herfindahl"))
  x2[, citypair:=NULL]
  x2 <- data.matrix(x2)
  
  gb2 <- sweep(x2, 1, v0, `*`)
  gb2 <- colMeans(gb2)
  
  gb <- c(gb1, gb2)
  
  gmm <- t(gb)%*%A%*%t(t(gb))
    
  x1a <<- x1
  x2a <<- x2
  pa <<- p
  v0a <<- v0 
  gba <<- gb
  
  return(gmm)
}