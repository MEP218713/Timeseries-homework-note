ga34<-matrix(c(0.43,-0.1),byrow = TRUE,nrow = 2)
ga2132<-matrix(c(0.23,-1.1,0.43,0.23),byrow = TRUE,ncol = 2)
ga2132
ga2132ni<-solve(ga2132)
a12<-ga2132%*%ga34
a012<-c(-1,0.2089,0.1619)
a012t<-matrix(c(-1,0.2089,0.1619),ncol=1)
gak0<-matrix(c(5.61,-1.1,0.23,
               -1.1,5.61,-1.1,
               0.23,-1.1,5.61),byrow=TRUE,ncol=3)
gak1<-matrix(c(-1.1,0.23,0.43,
               5.61,-1.1,0.23,
               -1.1,5.61,-1.1),byrow=TRUE,ncol=3)
gak2<-matrix(c(0.23,0.43,-0.1,
               -1.1,0.23,0.43,
               5.61,-1.1,0.23),byrow=TRUE,ncol=3)

gay0<-a012%*%gak0%*%a012t
gay1<-a012%*%gak1%*%a012t
gay2<-a012%*%gak2%*%a012t
gay0;gay1;gay2

## 从自协方差函数利用递推方法估计模型参数的函数
## Given \gamma_0, \gamma_1, \dots, \gamma_q,
## Solve MA(q) coefficients b_1, \dots, b_q, \sigma^2
## Using Li Lei's algorithm.
## Input: gms -- \gamma_0, \gamma_1, \dots, \gamma_q
ma.solve <- function(gms, k=100){
  q <- length(gms)-1
  if(q==1){
    rho1 <- gms[2] / gms[1]
    b <- (1 - sqrt(1 - 4*rho1^2))/(2*rho1)
    s2 <- gms[1] / (1 + b^2)
    return(list(b=b, s2=s2))
  }
  A <- matrix(0, nrow=q, ncol=q)
  for(j in seq(2,q)){
    A[j-1,j] <- 1
  }
  cc <- numeric(q); cc[1] <- 1
  gamma0 <- gms[1]
  gammas <- numeric(q+k)
  gammas[1:(q+1)] <- gms
  gamq <- gms[-1]
  Gammak <- matrix(0, nrow=k, ncol=k)
  for(ii in seq(k)){
    for(jj in seq(k)){
      Gammak[ii,jj] <- gammas[abs(ii-jj)+1]
    }
  }
  Omk <- matrix(0, nrow=q, ncol=k)
  for(ii in seq(q)){
    for(jj in seq(k)){
      Omk[ii,jj] <- gammas[ii+jj-1+1]
    }
  }
  PI <- Omk %*% solve(Gammak, t(Omk))
  s2 <- gamma0 - c(t(cc) %*% PI %*% cc)
  b <- 1/s2 * c(gamq - A %*% PI %*% cc)
  return(list(b=b, s2=s2))
}