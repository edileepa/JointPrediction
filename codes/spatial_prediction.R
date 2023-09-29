library(PrevMap)
library(rgdal)
library(geoR)
library(sp)
library(sf)
library(Matrix)
library(pdist)
library(splancs)
library(raster)


Rwanda_map <- st_read(".\\data\\RWA_adm\\RWA_adm0.shp")
#Rwanda_map <- st_read("D:\\OneDrive - Lancaster University\\github\\PredictionAlgorithm\\raw_data\\map.shp")

plot(st_geometry(Rwanda_map))
st_crs(Rwanda_map)
Rwanda_map.web <-st_transform(Rwanda_map,3857)


#import data
data <- read.csv(".\\data\\data.2008.cleaned.csv")
head(data)
dim(data)

ASCARIS.2008.estim     <- readRDS(".\\outputs\\ASCARIS.2008.estim.rds")
ASCARIS.2008.S.samples <- readRDS(".\\outputs\\ASCARIS.2008.S.samples.rds")
ls()

dim(data)
unique(data$ID)
dim(ASCARIS.2008.S.samples)

#step 2: spatial prediction---------

spatial_prediction_joint_model<-function(data,estim,S.samples,c1,c2,output ){
  
  p <- 1
  q <- 1
  
  beta1 <- estim$par[1:p]
  beta2 <- estim$par[(p+1):(p+q)] 
  
  sigma2.pi <- exp(estim$par[p+q+1]) 
  sigma2.lambda <-  exp(estim$par[p+q+2])
  
  phiS <-  exp(estim$par[p+q+3])
  phiT <-  exp(estim$par[p+q+4])
  
  gamma <- exp(estim$par[p+q+5])
  
  
  #ID.coords <- create.ID.coords(data,~LONGITUDE+LATITUDE) # ID of the unique set of coordinatesr
  
  data <- unique(data[,c("LATITUDE","LONGITUDE")])
  print(dim(data))
  
  coords.ll <- st_as_sf(data,coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(Rwanda_map))
  coords.web <-st_transform(coords.ll,3857)
  
  n.x <- nrow(st_coordinates(coords.web))
  cat("\nnumber of observed locations:",n.x,"\n")
  #n <- length(y)
  
  U <- dist(st_coordinates(coords.web)/1000) #distance in Km
  
  library(Matrix)
  cov.pars <- c(sigma2.pi,sigma2.lambda,phiS,phiT)
  cov.pars
  #cat("\nCovariance Par :",cov.pars,"\n")
  
  varcov <- function(cov.pars) {
    sigma2.pi     <- cov.pars[1]
    sigma2.lambda <- cov.pars[2]
    phiS          <- cov.pars[3]
    phiT          <- cov.pars[4]
    
    M <- matrix(NA,2*n.x,2*n.x)
    Sigma.shared <- varcov.spatial(dists.lowertri = U,
                                   cov.pars=c(1,phiT), #why sigma2 is 1? not sigma2.pi? reason is sigma2.pi is used below
                                   nugget = 0, kappa=0.5)$varcov
    M[1:n.x,1:n.x] <- sigma2.pi*(varcov.spatial(dists.lowertri = U,
                                                cov.pars=c(1,phiS), #why sigma2 is 1? not sigma2.lambda? reason is sigma2.lambda is used below
                                                nugget = 0, kappa=0.5)$varcov+Sigma.shared)
    
    M[(n.x+1):(2*n.x),(n.x+1):(2*n.x)] <- sigma2.lambda*(varcov.spatial(dists.lowertri = U,
                                                                        cov.pars=c(1,phiS),
                                                                        nugget = 0, kappa=0.5)$varcov+Sigma.shared)
    M[(n.x+1):(2*n.x),1:n.x] <- M[1:n.x,(n.x+1):(2*n.x)] <- sqrt(sigma2.pi*sigma2.lambda)*Sigma.shared
    M
  }
  
  #
  # library(geoR)
  M <- varcov(cov.pars)
  M1.inv <- solve(M[1:n.x,1:n.x])
  M2.inv <- solve(M[(n.x+1):(2*n.x),(n.x+1):(2*n.x)])
  
  # Create the prediction grid
  #poly <- coords[chull(coords),]
  
  #option 1:
  Rwanda_grid <- st_make_grid(Rwanda_map,cellsize = 0.1)
  Rwanda_grid <- Rwanda_grid[Rwanda_map] #remove this to get more area for prediction
  Rwanda_grid.web <- st_transform(Rwanda_grid,3857)
  Rwanda_grid.coords <- st_coordinates(Rwanda_grid.web)[,1:2]/1000

  #option 2:
  #Rwanda_grid.coords <- st_coordinates(Rwanda_map.web)[,1:2]/1000

  poly <- Rwanda_grid.coords[chull(Rwanda_grid.coords),]
  poly
  #grid.pred <- gridpts(poly,xs=5,ys=5)
  grid.pred <- gridpts(poly,npts=1880)

  
  par(mfrow=c(2,2))
  # plot(grid.pred)
  # plot(st_geometry(Rwanda_map.web)/1000,add=T)
  # points(st_coordinates(coords.web)/1000,col="red",pch=16)
  
  #library(pdist)
  #?pdist
  
  coords <- st_coordinates(coords.web)/1000
  U.pred.coords <- as.matrix(pdist(grid.pred,coords))
  dim(U.pred.coords)
  n.samples <- dim(S.samples)[1]
  
  # Predict S1+T
  C1 <- sigma2.pi*(exp(-U.pred.coords/phiS)+exp(-U.pred.coords/phiT))
  A1 <- as.matrix(C1%*%M1.inv)
  mean.cond.S1 <- sapply(1:n.samples, function(i) beta1+A1%*%S.samples[i,1:n.x]) # <-- this line
  
  #mu1 <- rep(beta1,nrow(S.samples))
  #mean.cond.S1 <- sapply(1:n.samples, function(i) beta1+A1%*%(S.samples[i,1:n.x]-mu1)) # <-- this does not give correct answer
  
  #sd.cond.S1 <- sqrt(sigma2.pi-apply(A1*C1,1,sum))
  sd.cond.S1 <- sqrt(2*sigma2.pi-apply(A1*C1,1,sum))
  S1.samples <- sapply(1:n.samples,function(i) mean.cond.S1[,i]+sd.cond.S1*rnorm(nrow(grid.pred)))
  
  S1.mean <- apply(S1.samples,1,mean)
  cat("\nsummary(S1.mean)",summary(S1.mean))
  
  # Predict S2+T
  C2 <- sigma2.lambda*(exp(-U.pred.coords/phiS)+exp(-U.pred.coords/phiT))
  A2 <- as.matrix(C2%*%M2.inv)
  mean.cond.S2 <- sapply(1:n.samples, function(i) beta2+A2%*%S.samples[i,(n.x+1):(2*n.x)]) # <-- this line
  
  #mu2 <- rep(beta2,nrow(S.samples))
  #mean.cond.S2 <- sapply(1:n.samples, function(i) beta2+A2%*%(S.samples[i,(n.x+1):(2*n.x)]-mu2)) # <-- this does not give correct answer
  
  #sd.cond.S2 <- sqrt(sigma2.lambda-apply(A2*C2,1,sum))
  sd.cond.S2 <- sqrt(2*sigma2.lambda-apply(A2*C2,1,sum))
  S2.samples <- sapply(1:n.samples,function(i) mean.cond.S2[,i]+sd.cond.S2*rnorm(nrow(grid.pred)))
  
  S2.mean <- apply(S2.samples,1,mean)
  cat("\nsummary(S2.mean)",summary(S2.mean))
  
  mean.pi.samples <- sapply(1:n.samples,
                            function(i) (exp(S1.samples[,i])/(1+exp(S1.samples[,i]))))
  
  mean.pi.mean <- apply(mean.pi.samples,1,mean)
  cat("\nsummary(mean.pi.mean)",summary(mean.pi.mean))
  
  
  mean.counts.samples <- sapply(1:n.samples,
                                function(i) (exp(S1.samples[,i])/(1+exp(S1.samples[,i])))*
                                  exp(S2.samples[,i])*gamma(1+1/gamma))
  
  mean.counts.mean <- apply(mean.counts.samples,1,mean)
  summary(mean.counts.mean) #minimum is 154 (but we want zeros)
  cat("\nsummary(mean.counts.mean)",summary(mean.counts.mean))
  #hist(mean.counts.mean)
  
  r.pi<- rasterFromXYZ(cbind(grid.pred,mean.pi.mean))
  plot(r.pi,main="Pi")
  points(coords,pch=20)
  plot(st_geometry(Rwanda_map.web)/1000,add=T)
  
  r<- rasterFromXYZ(cbind(grid.pred,mean.counts.mean))
  plot(r,main="Intensity")
  points(coords,pch=20)
  plot(st_geometry(Rwanda_map.web)/1000,add=T)
  
  mean.counts.prev.c1 <- apply(mean.counts.samples,1,function(x) mean(x>c1))
  r.prev.c1<- rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c1))
  plot(r.prev.c1,main=paste("Prev(Intensity >", c1, ")",sep=""))
  points(coords,pch=20)
  plot(st_geometry(Rwanda_map.web)/1000,add=T)
  
  mean.counts.prev.c2 <- apply(mean.counts.samples,1,function(x) mean(x>c2))
  r.prev.c2<- rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c2))
  plot(r.prev.c2,main=paste("Prev(Intensity >", c2, ")",sep=""))
  points(coords,pch=20)
  plot(st_geometry(Rwanda_map.web)/1000,add=T)
  
  #coordinate transfer to long and lat (otherwise cann't comapre with binom logistic model)
  # grid.pred.df <- data.frame(long=grid.pred[,1]*1000,
  #                            lat=grid.pred[,2]*1000,
  #                            pred=mean.counts.prev.c1)
  # 
  # grid.pred.sf <- st_as_sf(grid.pred.df,coords =c("long","lat"),crs =3857)
  # grid.pred.sf
  
  # par(mfrow=c(1,3))
  # plot(grid.pred)
  
  #grid.pred.sf <- st_transform(grid.pred.sf, crs = st_crs(Rwanda_map))
  # grid.pred.sf
  # plot(grid.pred.sf)
  # plot(st_coordinates(coords.ll))
  # plot(grid.pred.sf,add=T)
  # st_coordinates(grid.pred.sf)
  #plot(st_coordinates(grid.pred.sf))
  
  #dt <- cbind(st_coordinates(grid.pred.sf),grid.pred.sf$pred)
  #dt <- data.table(as.data.frame(dt, xy = TRUE))
  #setnames(dt, "p190001", "z")
  
  
  
  Rwanda_grid.coords <- st_coordinates(Rwanda_grid)[,1:2]
  poly <- Rwanda_grid.coords[chull(Rwanda_grid.coords),]
  grid.pred <- gridpts(poly,npts=1880)
  
  # plot(grid.pred)
  # plot(rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c1)),main="Output")
  # plot(st_geometry(Rwanda_map),add=T)
  
  par(mfrow=c(1,1))
  
  if(output=="None"){
    
    
  }else if(output=="Prevalence") {
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.pi.mean)),main="Output")
    plot(st_geometry(Rwanda_map),add=T)
    return(cbind(grid.pred,mean.pi.mean)) 
    
  } else if(output=="c1") {
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c1)),main="Output")
    plot(st_geometry(Rwanda_map),add=T)
    return(cbind(grid.pred,mean.counts.prev.c1))  
  }
  
  else if(output=="c2") {
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c2)),main="Output")
    plot(st_geometry(Rwanda_map),add=T)
    return(cbind(grid.pred,mean.counts.prev.c2))  
  }
  
  else if(output=="all") {
    
    par(mfrow=c(2,2))
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.counts.mean)),main="Intensity")
    points(coords,pch=20)
    plot(st_geometry(Rwanda_map.web)/1000,add=T)
    
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.pi.mean)),main="Prev(Intensity>0")
    plot(st_geometry(Rwanda_map),add=T)
    
    
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c1)),main="Prev(Intensity>5000")
    plot(st_geometry(Rwanda_map),add=T)
    
    plot(rasterFromXYZ(cbind(grid.pred,mean.counts.prev.c2)),main="Prev(Intensity>25000")
    plot(st_geometry(Rwanda_map),add=T)
    
    
  }
  
  
}

spatial_prediction_joint_model(data,
                               ASCARIS.2008.estim,
                               ASCARIS.2008.S.samples,
                               5000,50000,"None")
