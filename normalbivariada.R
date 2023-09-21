# Esta función gráfica la superficie de la FDP normal
# bivariada con sus curvas de nivel

library(plotly)
library(mvtnorm)

fdp_normal <- function(mu,sigma){
  densidad <- function(x1,x2){
    z <- dmvnorm(x=c(x1,x2),mean=mu,sigma=sigma,log=F)
    z
  }  
  densidad <- Vectorize(densidad)
  x <- seq(mu[1]-3*sqrt(sigma[1,1]),mu[1]+3*sqrt(sigma[1,1]),length.out=100)
  y <- seq(mu[2]-3*sqrt(sigma[2,2]),mu[2]+3*sqrt(sigma[2,2]),length.out=100)
  z <- outer(x,y,densidad)

df.list <- data.frame(x = seq(mu[1]-3*sqrt(sigma[1,1]),mu[1]+3*sqrt(sigma[1,1]),length.out=100),
                      y = seq(mu[2]-3*sqrt(sigma[2,2]),mu[2]+3*sqrt(sigma[2,2]),length.out=100),
                      z = outer(x,y,densidad))

fig <- plot_ly(df.list, x = x, y = y, z = z, type = "surface") %>% add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)
fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "Y1"),
    yaxis = list(title = "Y2"),
    zaxis = list(title = "FDP"),
    camera=list(
      eye = list(x=1.87, y=0.88, z=-0.64)
    )
  )
)

fig
}

if(FALSE){
  ####################
  ##### Ejemplos #####
  ####################
  
  #########################
  ## Covarianza negativa ##
  #########################
  mu1 <- c(3,7)
  sigma1 <- matrix(c(1.8,-0.6,-0.6,2.3),2,2)
  fdp_normal(mu=mu1,sigma=sigma1)
  #####################################
  ## Covarianza cero (independencia) ##
  #####################################
  mu2 <- c(3,7)
  sigma2 <- matrix(c(1.8,0,0,2.3),2,2)
  fdp_normal(mu=mu2,sigma=sigma2)
  #########################
  ## Covarianza positiva ##
  #########################
  mu3 <- c(-3,5)
  sigma3 <- matrix(c(4,0.6,0.6,2.3),2,2)
  fdp_normal(mu=mu3,sigma=sigma3)
}
