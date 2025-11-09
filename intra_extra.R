# Esta función permite establecer si un punto x0 es 
# de interpolación o extrapolación.

intra_extra <- function(x0,ajuste){
  x0 <- t(c(1,x0))
  X <- model.matrix(ajuste)
  h0 <- x0 %*% solve(t(X) %*% X) %*% t(x0)  
  hmax <- max(hatvalues(ajuste))
  if (h0 <= hmax) {
    cat("El punto es de interpolación (h0 =", round(h0, 5), ")\n")
  } else {
    cat("El punto es de extrapolación (h0 =", round(h0, 5), ")\n")
  }
}

if(FALSE){
  set.seed(5879)
  n <- 150
  x1 <- rnorm(n, mean = 50, sd = 10)   
  x2 <- rnorm(n, mean = 30, sd = 5)    
  x3 <- rnorm(n, mean = -10, sd = 8)
  x4 <- rnorm(n, mean = 0, sd = 6)
  error <- rnorm(n, mean = 0, sd = 8)
  y <- 5 + 0.8*x1 - 1.2*x2 + 2*x3 - x4 + error
  datos <- data.frame(y, x1, x2, x3, x4)
  ajuste <- lm(y ~ x1 + x2 + x3 + x4, data = datos)
  x0 <- c(12,30,-20,2)
  intra_extra(x0,ajuste)
  x0 <- c(40,25,-6,3)
  intra_extra(x0,ajuste)
}
