##This script is set to run remotely via a Linux bash script on a 36 core computer 
## Hence the "CommandArgs" lines

#args <- commandArgs(TRUE)
#Start <- as.numeric(args[1])
#Stop <- as.numeric(args[2])
#Start <- c(93, 95)


sim.Burrs <- dget(file = 'Thesis_Sim_10182018.txt')
for (i in Start) {
  source('RJMCMC.R')
  require(Rdistance)
  library(runjags)
  library(rjags)
All.Data <- data.frame(sim.Burrs[[i]]$all.burr)
Burrs <- data.frame(dist = All.Data$dist, size = (All.Data$diam)/10, Found = All.Data$found, 
                    Trans = All.Data$trans, Veg = All.Data$known_veg/100, Occ = All.Data$occ)
Found <- subset(Burrs, Burrs$Found == 1)
n.transects <- sim.Burrs[[i]]$n.lines
n <- nrow(Burrs)
ifelse(n == 200, nz <- 300, nz <- 450)


# require(Rdistance)
# detection.data <- data.frame(siteID = as.factor(Found$Trans), groupsize = 1, dist = Found$dist)
# transect.data <- site.data <-  data.frame(siteID = as.factor(c(1:n.transects)), length = sim.Burrs[[i]]$len)
# area = 10000
# est <- autoDistSamp(dist~1, detection.data, site.data, w.lo=0, w.hi=max(Found$dist),
#                     likelihoods = c("halfnorm", "uniform", "Gamma"),
#                     series = c("cosine", "hermite", "simple"), expansions =0:3, warn = TRUE,
#                     area = area, ci = 0.95, R = 500, plot.bs = FALSE, plot = FALSE)
# #est <- F.automated.CDA(detection.data, transect.data, w.lo=0, w.hi=max(Found$dist),
# #	likelihoods=c("halfnorm", "uniform", "Gamma"),
# #	series=c("cosine", "hermite", "simple"), expansions=0:3, warn=TRUE,
# #	area=area, ci=0.95, R=500, by.id=FALSE, plot.bs=FALSE, plot=FALSE)
# 
# pt.est.null <- est$n.hat
# ci.est.null <- unname(est$ci)
# truth <- (n/100 <= ci.est.null[2] & n/100 >= ci.est.null[1])

#######################
#bayesian method
modelstring.Vegetation = "
model
{
  for (i in 1:(nind +nz)) {
  w[i] ~ dbern(psi)					##augmentation 
  x[i] ~ dunif(0,Bx)					#distance from line for the missed ones; Bx = max(distances) 
  
  z[i] ~ dgamma(a[i],c.beta[i])T(4,65)		                
  a[i] <- shape[clust[i]]		## a depends on which cluster the size comes from
  c.beta[i] <- betaClust[clust[i]]	#precision is allowed to vary between clusters
  clust[i] ~ dcat( pClust[1:Nclust] )	## clusters are defined as categories, assuming Nclust clusters
  
  v[i] ~ dbeta(d,e)T(.05,.99)				## vegetation is unknown so beta prior
  sigma[i] <- exp(sigma.int+sigma.beta*z[i]+sigma.gamma*v[i])
  logp[i] <- -((x[i]*x[i])/(2*sigma[i]*sigma[i]))	#This is the normal distribution with an adjustment for covs
  p[i] <- exp(logp[i])*xi[i]
  xi[i] <- ifelse(z[i] < b.point, m*z[i]+intercept, 1) #if less than b.point, probability is linear; if larger than b.point, perfect detection on line
  
  mu[i] <- w[i]*p[i] 					## probabilty of seeing it for all the ones we DID see (where w[i] = 1)
  
  y[i] ~ dbern(mu[i])
  
  o[i]~ dbin(o2[i], 1)
  logit(o2[i]) <- o.int + z.beta*z[i]						
  }
  
  for (q in 1:q) {
  v2[q] ~ dbeta (d,e)T(.05,.99)		# all veg measurements at the site
  }
  
  p.online ~ dunif(0.4, 0.6)		#estimated detection on line for 5 cm burrows
  #b.point <- 20
  b.point ~ dunif(18, 22)
  m <- (1-p.online)/(b.point-5) 	#slope for detection on the line for smaller burrows		
  intercept <- p.online-(5*m)	## finding intercept via the detection of the 5cm burrow 
  
  
  sigma.int~ dnorm(0,s.int.tau)T(0,)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau)T(0,)
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau)T(,0)
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40)
  e~dunif(.1,40)
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) ## probability of each cluster = the probability of that category in the ddirch distribution
  
  
  psi~ dunif(0,1)			#exists or not		
  
  
  Occ <- sum(o)/(nind+nz)
  N <- sum(w)	
  D <- N/(2*L*Bx)
  Nt <- N*Occ
  Dt <- Nt/(2*L*Bx)	#tort density
  
  juvi1 <- sum(z < 13)/(nind+nz)
  juvi2 <- sum(z < 20)/(nind+nz)
  juvi3 <- sum(z < 10)/(nind+nz)
}
"

x <- c(Found$dist, rep(NA, nz))
y <- c(rep(1,nrow(Found)), rep(0, nz))
z <- c(Found$size, rep(NA, nz))
v <- c(Found$Veg, rep(NA, nz))
v <- ifelse(v < 0.01, .01, v)
v <- ifelse(v > 0.95, .95, v)
o <- c(Found$Occ, rep(NA, nz))
nind <- nrow(Found)

v2 <- sim.Burrs[[i]]$Vegcollect/100
v2 <- ifelse(v2 < 0.01, .01, v2)
v2 <- ifelse(v2 > 0.95, .95, v2)
q = length(v2)
Bx = max(Found$dist, na.rm = TRUE)
L <- sum(sim.Burrs[[i]]$len)*10^-4

# w.clust <- c(0.50, 0.1, 0.33, .03, .04)
# shape.clust <- c(10, 20, 30, 40, 66)
# rate.clust <- c(1,.7,1, 2, 2)
# y.clust <- Found$size
# Z <- do.call(cbind, lapply(1:5, function(j)
#   w.clust[j]*dgamma(y.clust, shape.clust[j], rate = rate.clust[j])))
# Z <- apply(Z, 1, function(x) which(x==max(x))[1])
# res <- mixgam.rjmcmc(y = y.clust, nsweep = 80000, kmax = 10, k = 5, w = w.clust, 
#                      shape = shape.clust, rate = rate.clust, Z,
#                      delta=1.5, xi=NULL, kappa=NULL, alpha=NULL,
#                      beta=NULL, g=NULL, h=NULL, verbose=TRUE)
# 
# ksave <- res$k.save
# groups <- round(table(ksave[-(1:40000)])/40000,4)
# Nclust <- as.numeric(names(groups)[which(groups == max(groups))])
# if(length(Nclust) != 1){Nclust <- sample(Nclust,1)}
# if(Nclust < 3){Nclust <- 3}

Nclust <- 4

clust = rep(NA,(nind +nz)) 	# no idea which cluster anything is in, so unknown	
clust[which.min(z)]=1 # smallest value assigned to cluster 1; generally represents juvis
clust[which.max(z)]=Nclust # highest value assigned to largest cluster value; generally represents large adults
w.i <- c(rep(1, nind), rep(0,nz))

### Run the model 
jd.test2 = list(nind= nind, nz = nz, L = L, Bx = Bx, x=x,y=y,z=z, Nclust = Nclust, q=q, o=o, v=v, 
                clust= clust, onesRepNclust = rep(1,Nclust))
ji.test2 <-list(w = w.i, betaClust =c(1.6, .8, 1.25, 1.2), pClust = c(.4, .3, .2, .1),
                sigma.beta = .025, sigma.int = 1, sigma.gamma = -1.5,
                shape = c(52, 9, 65, 55), z.beta = -0.016) 
jp.test2 <- c("D","Dt", "N","Nt", "betaClust", "pClust","shape","sigma.int", "sigma.beta", "p.online", 
              "b.point", "Occ", "sigma.gamma", "juvi1", "juvi2", "juvi3", "z.beta", "o.int")

Foo2 <- autorun.jags(modelstring.Vegetation, data = jd.test2, inits = list(ji.test2, ji.test2, ji.test2), monitor= jp.test2, 
                     startburnin = 10000, startsample = 40000, adapt = 1000, method = 'parallel', n.chain = 3, max.time="8h", 
                     psrf.target = 1.2,silent.jags = FALSE, summarise = TRUE)
Covmod2 <- summary(Foo2)



results = list()
results$simulation <- paste("Dataset #", i)
results$true <- c(paste('Truth'), n/(1000*1000*10^-4))
results$DISTANCE <- c(paste('Distance'), ci.est.null, pt.est.null)
results$FoundN <- nind
results$Mod <- paste('Veg_Model')
results$Covsum <- Covmod2

results
lapply(results, function(x) write.table(data.frame(x), 'LTDS_Simulation_Veg.csv', 
              append= T, sep=',', col.names = TRUE, row.names = TRUE ))
}
