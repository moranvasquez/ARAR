# Esta función permite calcular intervalos con una confianza
# del gamma*100% para cada coeficiente de regresión asociado
# al MRLN.

gamma <- 0.95
icreg <- function(ajuste,gamma){
  ajuste <- ajuste1
  coef <- unname(ajuste$coefficients)
  X <- model.matrix(ajuste)
  n <- nrow(X)
  r <- length(coef)-1
  alpha <- 1- gamma
  conf <- gamma*100
  const <- qt(1-alpha/2,n-r-1)
  a <- diag(solve(t(X)%*%X))
  MSE <- sum(ajuste1$residuals^2)/(n-r-1)
  nomfilas <- paste("beta_",0:r,sep="")
  infsup <- matrix(NA, nr=r+1,nc=2)
  for(j in 1:(r+1)){
    infsup[j,] <- c(coef[j] - const*sqrt(a[j]*MSE), coef[j] + const*sqrt(a[j]*MSE))
  }
  colnames(infsup) <- c("L.I","L.S")
  rownames(infsup) <- nomfilas
  cat("Intervalos del",conf,"%", "de confianza para los coeficientes de regresión","\n")
  return(infsup)
}

if(FALSE){
  library(readxl)
  antrop <- read_excel("body.xlsx")
  ajuste1 <- lm(weight~gender+thigh.girth, data=antrop)
  icreg(ajuste1,0.95)
  icreg(ajuste1,0.99)
}
