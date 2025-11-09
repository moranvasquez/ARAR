# Esta función permite graficar la región (elipse) de las variables
# explicativas asociado a un MRLN con dos variables explicativas. Además
# ubica el punto x0=(x1,x2)' para determinar si es un punto de interpolación
# o extrapolación.

intra_extra_plot <- function(x1,x2,ajuste){
  if (!inherits(ajuste, "lm")) stop("'Ajuste' debe ser un objeto 'lm'.")
  X <- model.matrix(ajuste)
  if (ncol(X) != 3) stop("La región se grafica para el caso de dos variables explicativas.")
  S <- solve(t(X) %*% X)  # S = (X'X)^{-1}
  s11  <- S[1,1]
  s12  <- matrix(S[1, -1], ncol = 1)   
  S22  <- S[-1, -1]                    
  hmax <- max(hatvalues(ajuste))
  mu <- - solve(S22, s12)   # 2x1
  denom <- as.numeric(hmax - s11 + t(s12) %*% solve(S22, s12))
  A <- S22 / denom          # 2x2
  eig <- eigen(A)
  vals <- eig$values
  vecs <- eig$vectors
  radios <- 1 / sqrt(vals)  
  theta <- seq(0, 2*pi, length.out = 100)
  unit_circle <- cbind(cos(theta), sin(theta))   # npoints x 2
  ellipse_pts <- t( vecs %*% (diag(radios) %*% t(unit_circle)) )  # npoints x 2
  ellipse_pts <- sweep(ellipse_pts, 2, as.numeric(mu), FUN = "+")
  datos <- as.data.frame(X[, -1])
  names(datos) <- colnames(X)[-1]
  
  x.lim <- c(mean(datos[,1]) - 6*sd(datos[,1]),mean(datos[,1]) + 6*sd(datos[,1]))
  y.lim <- c(mean(datos[,2]) - 6*sd(datos[,2]),mean(datos[,2]) + 6*sd(datos[,2]))
  plot(datos[,1], datos[,2], pch = 19,
       xlab = expression(x[1]), ylab = expression(x[2]),
       xlim = x.lim, ylim = y.lim, col="gray",
       main=bquote("Región de las variables explicativas:  " ~ h[max] == .(sprintf("%.4f", hmax))))
  lines(ellipse_pts[,1], ellipse_pts[,2], col = "blue", lwd = 2)
  points(x1, x2, pch = 19, col = "red")
  
  x0 <- t(c(1,x1,x2))
  h0 <- x0 %*% solve(t(X) %*% X) %*% t(x0)
  if (h0 <= hmax) {
    cat("El punto es de interpolación (h0 =", round(h0, 5), ")\n")
  } else {
    cat("El punto es de extrapolación (h0 =", round(h0, 5), ")\n")
  }
}

if(FALSE){
  library(mnormt)
  set.seed(5879)
  n <- 150
  mu <- c(50,30)
  S <- matrix(c(3,-2,-2,2),2,2)
  Xr <- rmnorm(n = n, mean = mu, varcov=S)
  error <- rnorm(n, mean = 0, sd = 8)
  y <- 5 + 0.8*Xr[,1] - 1.2*Xr[,2] + error
  datos <- data.frame(y, Xr[,1], Xr[,2])
  names(datos) <- c("y","x1","x2")
  ajuste <- lm(y ~ x1 + x2, data = datos)
  intra_extra_plot(47,25,ajuste)
  intra_extra_plot(52,28,ajuste)
}

