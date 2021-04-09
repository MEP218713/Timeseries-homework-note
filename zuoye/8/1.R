install.packages("BB")
library(BB)

fun <- function(x) 
  { 
  f <- numeric(length(x)) 					
  f[1] <-  x[3]*(x[1]*x[1]+x[2]*x[2])-12.4168
  f[2] <-  x[3]*(x[1]*x[2]+x[1])+4.7520
  f[3] <-  x[3]*x[2]-5.2
  f 
} 
startx <- c(1,2,3)
result = dfsane(startx,fun,control=list(maxit=2500,trace = FALSE))
theta = result$par