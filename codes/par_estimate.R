
library(geoR)
library(PrevMap)
library(sf)


#step 1: parameter estimate + generate S.sample ------

parameter.estimate <- function(data,y,estim,hess){
  
  ID.coords <- create.ID.coords(data,~LONGITUDE+LATITUDE) # ID of the unique set of coordinatesr
  
  par(mfrow=c(1,2))
  # coords.ll <- SpatialPoints(unique(data[,c("LONGITUDE","LATITUDE")]),CRS("+init=epsg:4236")) #coordinates in longitude and latitude
  # plot(coordinates(coords.ll),main="Coordinates in degrees") #coordinates in degrees
  # coords.web <- spTransform(coords.ll,CRS("+init=epsg:3857")) #coordinates in web (in meters)
  # plot(coordinates(coords.web),main="Coordinates in meters")
  
  data <- unique(data[,c("LATITUDE","LONGITUDE")])
  coords.ll <- st_as_sf(data,coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(Rwanda_map))
  coords.web <-st_transform(coords.ll,3857)
  
  
  p <- 1
  q <- 1
  
  beta1.0 <- estim$par[1:p]         # mean of linear predictors of prevalence (binomial logistic model - intercept only)
  beta2.0 <- estim$par[(p+1):(p+q)] # mean of linear predictors of intensity (log-linear model - intercept only)
  
  sigma2.pi0     <- exp(estim$par[p+q+1])  # variance of linear predictors of prevalence (binomial logistic model)
  sigma2.lambda0 <-  exp(estim$par[p+q+2]) # variance of linear predictors of intensity (log-linear model)
  
  phiS.0 <-  exp(estim$par[p+q+3]) # scale of the spatial correlation --> range of the spatial correlation
  phiT.0 <-  exp(estim$par[p+q+4]) # scale of the spatial correlation --> range of the spatial correlation
  
  gamma0 <- exp(estim$par[p+q+5]) # Shape parameter Weibull
  c(beta1.0,beta2.0,sigma2.pi0 ,sigma2.lambda0,phiS.0,phiT.0,gamma0)
  
  #y <- data$mfcount
  #y <- data.2008$anySTH
  
  ind.pos <- which(y>0)
  n.zero.cases <- tapply(y,ID.coords,function(x) length(which(x==0)))
  n.pos.cases <- tapply(y,ID.coords,function(x) length(which(x>0)))
  
  n.x <- nrow(st_coordinates(coords.web))
  n <- length(y)
  U <- dist(st_coordinates(coords.web)/1000)
  
  cov.pars0 <- c(sigma2.pi0,sigma2.lambda0,phiS.0,phiT.0)
  cov.pars0
  
  
  varcov <- function(cov.pars) {
    sigma2.pi <- cov.pars[1]
    sigma2.lambda <- cov.pars[2]
    phiS <- cov.pars[3]
    phiT <- cov.pars[4]
    
    M <- matrix(NA,2*n.x,2*n.x)
    Sigma.shared <- varcov.spatial(dists.lowertri = U,
                                   cov.pars=c(1,phiT),
                                   nugget = 0, kappa=0.5)$varcov
    M[1:n.x,1:n.x] <- sigma2.pi*(varcov.spatial(dists.lowertri = U,
                                                cov.pars=c(1,phiS),
                                                nugget = 0, kappa=0.5)$varcov+Sigma.shared)
    M[(n.x+1):(2*n.x),(n.x+1):(2*n.x)] <- sigma2.lambda*(varcov.spatial(dists.lowertri = U,
                                                                        cov.pars=c(1,phiS),
                                                                        nugget = 0, kappa=0.5)$varcov+Sigma.shared)
    M[(n.x+1):(2*n.x),1:n.x] <- M[1:n.x,(n.x+1):(2*n.x)] <- sqrt(sigma2.pi*sigma2.lambda)*Sigma.shared
    M
  }
  
  M0 <- varcov(cov.pars0)
  M0.inv <- solve(M0) #get the inverse
  
  ind.S2 <- ID.coords[ind.pos]
  ind.loc.pos <- which(tapply(y,ID.coords,sum)>0) #cases more than zero
  
  integrand <- function(S.tot) {
    S1 <- S.tot[1:n.x]
    S2 <- S.tot[(n.x+1):(2*n.x)]
    
    rand.eff.llik <- -0.5*t(S.tot)%*%M0.inv%*%S.tot
    eta.prob <- beta1.0+S1
    eta.counts <- beta2.0+S2[ind.S2]
    prob <- exp(eta.prob)/(1+exp(eta.prob))
    lambda <- exp(eta.counts)
    llik.weib <- sum(n.zero.cases*log(1-prob))+
      sum(n.pos.cases*log(prob))+
      sum((gamma0-1)*log(y[ind.pos]/lambda)+
            -(y[ind.pos]/lambda)^gamma0+
            log(gamma0/lambda))
    
    out <- rand.eff.llik+llik.weib
    as.numeric(out)
  }
  
  grad.integrand <- function(S.tot) {
    S1 <- S.tot[1:n.x]
    S2 <- S.tot[(n.x+1):(2*n.x)]
    
    eta.prob <- beta1.0+S1
    eta.counts <- beta2.0+S2[ind.S2]
    
    g <- as.numeric(-M0.inv%*%S.tot)
    
    prob <- exp(eta.prob)/(1+exp(eta.prob))
    der.prob <- prob/(1+exp(eta.prob))
    lambda <- exp(eta.counts)
    
    g[1:n.x] <- g[1:n.x]-(n.zero.cases/(1-prob))*der.prob
    g[1:n.x] <- g[1:n.x]+n.pos.cases*der.prob/prob
    
    g[(n.x+1):(2*n.x)][ind.loc.pos] <- g[(n.x+1):(2*n.x)][ind.loc.pos]+
      tapply(-gamma0-(gamma0*(y[ind.pos]/lambda)^(gamma0-1))*
               (-y[ind.pos]/(lambda^2))*lambda,ind.S2,sum)
    as.numeric(g)
  }
  
  
  maxim.integrand <- maxBFGS(
    integrand,
    grad.integrand,
    start=rep(0,2*n.x),         
    print.level = 1)
  
  
  Sigma.tilde <- solve(-maxim.integrand$hessian)
  Sigma.tilde.sroot <- t(chol(Sigma.tilde))
  A <- solve(Sigma.tilde.sroot)
  
  n.sim <- 110000
  burnin <- 10000
  thin <- 10
  
  S.hat <- maxim.integrand$estimate
  Sigma.W <- A%*%M0%*%t(A)
  Sigma.W.inv <- solve(Sigma.W) #get the inverse
  mean.W <- as.numeric(-A%*%S.hat)
  
  
  log.dens <- function(W,S.tot) {
    S1 <- S.tot[1:n.x]
    S2 <- S.tot[(n.x+1):(2*n.x)]
    
    rand.eff.llik <- -0.5*t(W-mean.W)%*%Sigma.W.inv%*%(W-mean.W)
    eta.prob <- beta1.0+S1
    eta.counts <- beta2.0+S2[ind.S2]
    prob <- exp(eta.prob)/(1+exp(eta.prob))
    lambda <- exp(eta.counts)
    llik.weib <- sum(n.zero.cases*log(1-prob))+
      sum(n.pos.cases*log(prob))+
      sum((gamma0-1)*log(y[ind.pos]/lambda)+
            -(y[ind.pos]/lambda)^gamma0+
            log(gamma0/lambda))
    
    out <- rand.eff.llik+llik.weib
    as.numeric(out)
  }
  
  grad.log.dens <- function(W,S.tot) {
    S1 <- S.tot[1:n.x]
    S2 <- S.tot[(n.x+1):(2*n.x)]
    
    eta.prob <- beta1.0+S1
    eta.counts <- beta2.0+S2[ind.S2]
    
    g.W <- as.numeric(-Sigma.W.inv%*%(W-mean.W))
    
    prob <- exp(eta.prob)/(1+exp(eta.prob))
    der.prob <- prob/(1+exp(eta.prob))
    lambda <- exp(eta.counts)
    
    g <- rep(0,2*n.x)
    g[1:n.x] <- g[1:n.x]-(n.zero.cases/(1-prob))*der.prob
    g[1:n.x] <- g[1:n.x]+n.pos.cases*der.prob/prob
    
    g[(n.x+1):(2*n.x)][ind.loc.pos] <- g[(n.x+1):(2*n.x)][ind.loc.pos]+
      tapply(-gamma0-(gamma0*(y[ind.pos]/lambda)^(gamma0-1))*
               (-y[ind.pos]/(lambda^2))*lambda,ind.S2,sum)
    as.numeric(g.W+t(Sigma.tilde.sroot)%*%g) 
  }
  
  
  W.curr <- rep(0,2*n.x)
  S.tot.curr <- as.numeric(Sigma.tilde.sroot%*%W.curr+S.hat)
  lp.curr <- log.dens(W.curr,S.tot.curr)
  h <- 1.65/((2*n.x)^(1/6))
  h.vec <- rep(NA,n.sim)
  acc <- 0
  
  mean.curr <- as.numeric(W.curr + (h^2/2)*grad.log.dens(W.curr,S.tot.curr))
  
  n.samples <- (n.sim-burnin)/thin
  S.samples <- matrix(NA,n.samples,2*n.x)
  
  for(i in 1:n.sim) {
    W.prop <- mean.curr+h*rnorm(2*n.x)
    S.tot.prop <-  as.numeric(Sigma.tilde.sroot%*%W.prop+S.hat)
    mean.prop <- as.numeric(W.prop + (h^2/2)*grad.log.dens(W.prop,S.tot.prop))
    lp.prop <- log.dens(W.prop,S.tot.prop)
    
    dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2))
    dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
    
    log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr
    
    if(log(runif(1)) < log.prob) {
      acc <- acc+1
      W.curr <- W.prop
      S.tot.curr <- S.tot.prop
      lp.curr <- lp.prop
      mean.curr <- mean.prop
    }
    
    if( i > burnin & (i-burnin)%%thin==0) {
      S.samples[(i-burnin)/thin,] <- S.tot.curr
    }
    
    h.vec[i] <- h <- max(0,h + 0.01*i^(-0.001)*(acc/i-0.57))
    cat("iteration:",i,"\r")
  } 
  
  
  par0 <- c(beta1.0,beta2.0,log(c(sigma2.pi0,sigma2.lambda0,phiS.0,phiT.0)),log(gamma0))
  #par0
  
  
  
  compute.llik <- function(S.tot,par,log.det.M,M.inv){
    S1 <- S.tot[1:n.x]
    S2 <- S.tot[(n.x+1):(2*n.x)]
    
    p <- 1
    q <- 1
    beta1 <- par[1]  #these need to be changed 
    beta2 <- par[2]  #these need to be changed 
    gamma <- exp(par[p+q+5])
    rand.eff.llik <- -0.5*(log.det.M+t(S.tot)%*%M.inv%*%S.tot)
    eta.prob <- beta1+S1              #eta.prob <- D1%*%beta1+S1 for covariates
    eta.counts <- beta2+S2[ind.S2]    #eta.counts <- D2%*%beta2+S2 for covariates
    prob <- exp(eta.prob)/(1+exp(eta.prob))
    lambda <- exp(eta.counts)
    llik.weib <- sum(n.zero.cases*log(1-prob))+
      sum(n.pos.cases*log(prob))+
      sum((gamma-1)*log(y[ind.pos]/lambda)+
            -(y[ind.pos]/lambda)^gamma+
            log(gamma/lambda))
    
    out <- rand.eff.llik+llik.weib
    as.numeric(out)
  }
  
  p <- 1
  q <- 1
  D1 <- cbind(rep(1,nrow(data)))
  D2 <- cbind(rep(1,nrow(data)))
  
  log.MC.lik <- function(par) {
    sigma2.pi <- exp(par[p+q+1])
    sigma2.lambda <- exp(par[p+q+2])
    phiS <- exp(par[p+q+3])
    phiT <- exp(par[p+q+4])
    cov.pars <- c(sigma2.pi,sigma2.lambda,phiS,phiT)
    M <- varcov(cov.pars)
    M.inv <- solve(M)
    log.det.M <- as.numeric(determinant(M)$modulus)
    res <- sapply(1:n.samples,function(i) compute.llik(S.samples[i,],par,log.det.M,M.inv))
  }
  
  log.f.tilde <- log.MC.lik(par0)
  length(log.f.tilde)
  log.f.tilde[1:10]
  
  objective.f <- function(par) {
    log.f <- log.MC.lik(par)
    log(mean(exp(log.f-log.f.tilde)))
  }
  
  
  
  estim <- nlminb(par0,
                  function(x) -objective.f(x),
                  control=list(trace=1))
  estim
  cat("\nParameter estimates:", estim$par)
  
  if(hess =="yes"){
    
    library(numDeriv)
    cat("\nCalculating hessian...........")
    estim$hessian <- -hessian(objective.f,estim$par)
    estim$hessian
    save(estim,file=paste(name, "geo_weibull.RData",sep="_"))
    
    cat("\nCalculating gradient...........")
    estim$gradient <- -grad(objective.f,estim$par)
    estim$gradient
    save(estim,file=paste(name, "geo_weibull.RData",sep="_"))
    
    
  }
  
  #estim
  
  #save(S.samples,file=paste(name,"_geo_weibull_samples.RData",sep="")) #this was missing in the original code
  #cat("\n S.Samples saved under the name.....",paste(name,"_geo_weibull_samples.RData",sep=""),"\n")
  
  return(list(estim,S.samples))
  
  
}

#assign this output
parameter.estimate(data,y,estim,hess)



head(data)
ASCARIS.estim <- data.frame(par=c(-1.2273226,  8.5273829,  3.5659694,  1.5593032 , 5.6983150,  6.4016951, -0.1888495))
ASCARIS.estim 

for(i in 1:3){
  print(ASCARIS.estim$par)
  store <- parameter.estimate(data =data,
                              y    =data$ASCARISAVERAGE,
                              estim=ASCARIS.estim,
                              hess="no")
  ASCARIS.estim <- store[[1]]
  ASCARIS.S.samples <- store[[2]]
  print(str(ASCARIS.estim) )
  print(str(ASCARIS.S.samples) )
  print(ASCARIS.estim$par)
  print(ASCARIS.estim$objective)
  print(ASCARIS.estim)
}

save.image()
