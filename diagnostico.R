# Esta función elabora gráficos de residuales
# para evaluar la calidad de ajuste del MRLN

#=============================================
# Guargar el MRLN en un objeto llamado ajuste
#=============================================

#=========
# Outliers
#=========
resid_ti <- rstudent(ajuste)
eje.x <- round(quantile(1:n, probs = c(0, 0.5, 1)), 0)
eje.y <- c(-3,0,3)
ylim <- range(-3,3,resid_ti)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
plot(resid_ti,xaxt="n",yaxt="n",xlab="",ylab="",pch=19,ylim=ylim)
abline(h=0,col="blue",lty = "dashed")
abline(h=3,col="red",lty = "dashed")
abline(h=-3,col="red",lty = "dashed")
axis(side = 1, at = eje.x, cex.axis=1.8,line = 1,tick = FALSE)
axis(side = 2, at = eje.y, cex.axis=1.8,line = -0.4,tick = FALSE)
title(xlab="Índice", line=7, cex.lab=2, cex=1.5)
title(ylab="Residual estud. exter.", line=5, cex.lab=2, cex=1.5)
posibles_outliers <- unname(which(abs(resid_ti) > 3))
if (length(posibles_outliers) == 0){
  cat("No se detectaron posibles outliers (|t_i| > 3).\n")
} else{
  resultados <- data.frame(posibles_outliers,resid_ti[posibles_outliers])
  rownames(resultados) <- NULL
  names(resultados) <- c("Observación", "t_i")
  cat("Posibles outliers detectados (|t_i| > 3):\n")
  print(resultados, row.names = FALSE)
}

#==================
# Puntos de palanca
#==================
X <- model.matrix(ajuste)
n <- nrow(X)
r <- ncol(X)-1
hii <- hatvalues(ajuste)
hmax <- max(hii)
lim <- 2*(r+1)/n
eje.x <- round(quantile(1:n, probs = c(0, 0.5, 1)), 0)
eje.y <- round(seq(min(hii), max(hii), length.out = 3), 3)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
plot(hii,xaxt="n",yaxt="n",xlab="",ylab="",pch=19)
axis(side = 1, at = eje.x, cex.axis=1.8,line = 1,tick = FALSE)
axis(side = 2, at = eje.y, cex.axis=1.8,line = -0.4,tick = FALSE)
title(xlab="Índice", line=7, cex.lab=2, cex=1.5)
title(ylab="Leverage", line=5, cex.lab=2, cex=1.5)
abline(h=lim,col="red",lty = "dashed")
posibles_palanca <- unname(which(abs(hii) > lim))
if (length(posibles_palanca) == 0){
  cat(paste0("No se detectaron posibles puntos de palanca (|h_ii| > ", lim, ").\n"))
} else{
  resultados <- data.frame(posibles_palanca,hii[posibles_palanca])
  rownames(resultados) <- NULL
  names(resultados) <- c("Observación", "h_ii")
  cat(paste0("Posibles puntos de palanca detectados (|h_ii| > ", lim, ").\n"))
  print(resultados, row.names = FALSE)
}

#===================
# Puntos influyentes
#===================
X <- model.matrix(ajuste)
n <- nrow(X)
r <- ncol(X)-1
Di <- cooks.distance(ajuste)
lim.f <- qf(0.5,r+1,n-r-1)
lim.max <- max(Di,lim.f)
eje.x <- round(quantile(1:n, probs = c(0, 0.5, 1)), 0)
eje.y <- round(seq(min(resid_ti), max(resid_ti), length.out = 3), 1)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
plot(Di,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,lim.max),pch=19)
axis(side = 1, at = eje.x, cex.axis=1.8,line = 1,tick = FALSE)
axis(side = 2, at = eje.y, cex.axis=1.8,line = -0.4,tick = FALSE)
title(xlab="Índice", line=7, cex.lab=2, cex=1.5)
title(ylab="Distancia de Cook", line=5, cex.lab=2, cex=1.5)
abline(h=lim.f,col="red",lty = "dashed")
posibles_influyentes <- unname(which(Di > lim.f))
if (length(posibles_influyentes) == 0){
  cat(paste0("No se detectaron posibles puntos de influyentes (D_i > ", lim.f, ").\n"))
} else{
  resultados <- data.frame(posibles_influyentes,Di[posibles_influyentes])
  rownames(resultados) <- NULL
  names(resultados) <- c("Observación", "D_i")
  cat(paste0("Posibles puntos de influyentes detectados (|D_i| > ", lim, ").\n"))
  print(resultados, row.names = FALSE)
}

#================================
# Heterocedasticidad y linealidad
#================================
ajus_i <- ajuste$fitted.values
resid_ti <- rstudent(ajuste)
eje.x <- round(seq(min(ajus_i), max(ajus_i), length.out = 4), 1)
eje.y <- round(seq(min(resid_ti), max(resid_ti), length.out = 4), 1)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
plot(ajus_i,resid_ti,xaxt="n",yaxt="n",xlab="",ylab="",pch=19)
axis(side = 1, at = eje.x, cex.axis=1.8,line = 1,tick = FALSE)
axis(side = 2, at = eje.y, cex.axis=1.8,line = -0.4,tick = FALSE)
title(xlab="Valor ajustado", line=7, cex.lab=2, cex=1.5)
title(ylab="Residual estud. exter.", line=5, cex.lab=2, cex=1.5)

#============
# Normalidad
#============
X <- model.matrix(ajuste)
n <- nrow(X)
r <- ncol(X)-1
resid_ti <- sort(rstudent(ajuste))
ident =   diag(n)
H <- X %*% solve(t(X) %*% X) %*% t(X)
epsilon =   matrix(0,n,100)
e =   matrix(0,n,100)
e1 =   numeric(n)
e2 =   numeric(n)
for(i in 1:100){
  epsilon[,i] =   rnorm(n,0,1)
  e[,i] =   (ident - H)%*%epsilon[,i]
  u =   diag(ident - H)
  e[,i] =   e[,i]/sqrt(u)
  e[,i] =   sort(e[,i]) }
for(i in 1:n){
  eo =   sort(e[i,])
  e1[i] =   min(eo)
  e2[i] =   max(eo) }
medias <- rowMeans(e)
par(mfrow=c(1,1),pty="s",mex=0.4)
par(mar=c(10,10,10,5))
probs <- 1:n/(n+1)
cuants <- sort(qt(probs,df=n-r-2))
xlim = c(min(cuants), max(cuants))
ylim = range(e1,e2,resid_ti)
eje.x <- round(seq(min(cuants), max(cuants), length.out = 4), 1)
eje.y <- round(seq(min(resid_ti), max(resid_ti), length.out = 4), 1)
plot(cuants,resid_ti,xaxt="n",yaxt="n",xlab="",ylab="",pch=19,cex=0.7,xlim=xlim,ylim=ylim)
par(new=TRUE)
plot(cuants,e1,xaxt="n",yaxt="n",xlab="",ylab="",type="l",col="red",xlim=xlim,ylim=ylim)
par(new=TRUE)
plot(cuants,e2,xaxt="n",yaxt="n",xlab="",ylab="",type="l",col="red",xlim=xlim,ylim=ylim)
par(new=TRUE)
plot(cuants,medias,xaxt="n",yaxt="n",xlab="",ylab="",type="l",lty=2,col="blue",xlim=xlim,ylim=ylim)
axis(side = 1, at = eje.x, cex.axis=1.8,line = 1,tick = FALSE)
axis(side = 2, at = eje.y, cex.axis=1.8,line = -0.4,tick = FALSE)
title(xlab="Cuantil teórico", line=7, cex.lab=2, cex=1.5)
title(ylab="Cuantil observado", line=5, cex.lab=2, cex=1.5)
puntos.fuera <- sum(sort(resid_ti) < sort(e1) | sort(resid_ti) > sort(e2))
porc.puntos.fuera <- round(puntos.fuera*100/n,2)
cat(paste0(porc.puntos.fuera, "% de los puntos se salen de la envolvente simulada\n"))
