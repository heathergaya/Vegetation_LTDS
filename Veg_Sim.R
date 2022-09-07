##This script can be set to run remotely via a Linux bash script on a 36 core computer 
## Hence the "CommandArgs" lines
#Use the following if you want to run these through Linux:
#args <- commandArgs(TRUE) 
#Start <- as.numeric(args[1])
#Stop <- as.numeric(args[2])


sim.Burrs <- dget(file = 'Thesis_Sim_10182018.txt') #simulated data set
Start <- 1 #the first simulation to run
Stop <- nrow(summary(sim.Burrs)) #the last simulation you want to run
for (i in Start:Stop) {  #loop across each simulation
  require(Rdistance) #for comparing w/ maximum likelihood estimation
  library(runjags)
All.Data <- data.frame(sim.Burrs[[i]]$all.burr) #grab the simulation data 
Burrs <- data.frame(dist = All.Data$dist, size = (All.Data$diam)/10, Found = All.Data$found, 
                    Trans = All.Data$trans, Veg = All.Data$known_veg/100, Occ = All.Data$occ)
#need distance to transect, size of burrow, if the burrow was detected, what transect (for Rdistance), the veg information at that burrow and the occupancy of the burrow 

Found <- subset(Burrs, Burrs$Found == 1) #we only work with detected burrows 
n.transects <- sim.Burrs[[i]]$n.lines
n <- nrow(Burrs)
ifelse(n == 200, nz <- 300, nz <- 450) #data augmentation 

## If you want to compare with Program DISTANCE:
# require(Rdistance)
# detection.data <- data.frame(siteID = as.factor(Found$Trans), groupsize = 1, dist = Found$dist)
# transect.data <- site.data <-  data.frame(siteID = as.factor(c(1:n.transects)), length = sim.Burrs[[i]]$len)
# area <- 10000
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
  
  z[i] ~ dgamma(a[i],c.beta[i])T(5,65)	#burrow sizes	                
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
  logit(o2[i]) <- o.int + z.beta*z[i]			#larger burrows probably less likely to be occupied		
  ow[i] <- o[i]*w[i] # was this burrow both occupied and real
  }
  
  for (q in 1:q) {
  v2[q] ~ dbeta (d,e)T(.05,.99)		# all veg measurements at the site inform the veg at augmented burrows
  }
  
  p.online ~ dunif(.2, 0.8)		#estimated detection on line for 5 cm burrows
  #b.point <- 20 #could alternatively set break-point for perfect detection 
  b.point ~ dunif(15, 25) #prior for size of burrow where perfect detection is reached
  m <- (1-p.online)/(b.point-5) 	#slope for detection on the line for smaller burrows		
  intercept <- p.online-(5*m)	## finding intercept via the detection of the 5cm burrow 
  
  ## Priors are included below: 
  
  sigma.int~ dnorm(0,s.int.tau) #intercept of sigma (the shape parameter detection)
  s.int.tau <- 1/(s.int.sd*s.int.sd)
  s.int.sd ~ dunif(.00001,5)	
  sigma.beta~ dnorm(0,s.beta.tau) #relationship of sigma to burrow size
  s.beta.tau <- 1/(s.beta.sd*s.beta.sd)
  s.beta.sd ~ dunif(.00001,5)
  o.int~ dnorm(0,o.int.tau) #intercept for burrow occupancy (logit scale)
  o.int.tau <- 1/(o.int.sd*o.int.sd)
  o.int.sd ~ dunif(.00001,5)		
  z.beta~ dnorm(0,z.beta.tau) #relationship burrow occupancy to size (logit scale)
  z.beta.tau <- 1/(z.beta.sd*z.beta.sd)
  z.beta.sd ~ dunif(.00001,5)	
  sigma.gamma~ dnorm(0,s.gamma.tau) #relationship between vegetation and sigma parameter of detection 
  s.gamma.tau <- 1/(s.gamma.sd*s.gamma.sd)
  s.gamma.sd ~ dunif(.00001,5)
  
  d~dunif(.1,40) #for veg
  e~dunif(.1,40) # for veg 
  
  for (clustIdx in 1: Nclust) {
  shape[clustIdx] ~ dunif(1,80)
  betaClust[clustIdx] ~ dunif(0.2,2)
  }
  
  pClust[1:Nclust] ~ ddirch(onesRepNclust) ## probability of each cluster = the probability of that category in the ddirch distribution
  
  psi~ dunif(0,1)			#exists or not		
  
  # Derived Quantities
  
  Occ <- sum(ow[])/(nind+nz) #proportion occupied of all indiviudals (real and not)
  N <- sum(w)	 #Burrow Abundance 
  D <- N/(2*L*Bx) #density of burrows 
  Nt <- sum(ow[]) #occupied burrow abundance 
  Dt <- Nt/(2*L*Bx)	#tortoise (occupied burrow) density
  
  juvi1 <- sum(ow*z < 13)/Nt   #proportion burrows < 13 cm wide
  juvi2 <-  (sum(ow*z < 21)- sum(wo*z <13))/Nt #proportion burrows between 13 and 21 cm wide
  juvi3 <- sum(ow*z >= 21)/Nt #proportion burrows > 21 cm wide
}
"

x <- c(Found$dist, rep(NA, nz)) #distance to transect
y <- c(rep(1,nrow(Found)), rep(0, nz)) #found (1) or not found (0)
z <- c(Found$size, rep(NA, nz)) # size of burrows
v <- c(Found$Veg, rep(NA, nz)) #vegetation obstruction at burrows
v <- ifelse(v < 0.01, .01, v) #minimize extremes to help with convergence
v <- ifelse(v > 0.95, .95, v) #minimize extremes to help with convergence
o <- c(Found$Occ, rep(NA, nz)) #occupancy of burrows 
nind <- nrow(Found) #how many burrows were observed

v2 <- sim.Burrs[[i]]$Vegcollect/100 #vegetation everywhere at the site, not just at burrows 
v2 <- ifelse(v2 < 0.01, .01, v2) #minimize extremes to help with convergence
v2 <- ifelse(v2 > 0.95, .95, v2) #minimize extremes to help with convergence
q <- length(v2) #how many veg points were taken
Bx <- max(Found$dist, na.rm = TRUE) #maximum distance from transect to consider in abundance ests
L <- sum(sim.Burrs[[i]]$len)*10^-4 #length of transects

Nclust <- 5 #a large number of clusters

clust <- rep(NA,(nind +nz)) 	# no idea which cluster anything is in, so unknown	
clust[which.min(z)] <- 1 # smallest value assigned to cluster 1; generally represents juvis
clust[which.max(z)] <- Nclust # highest value assigned to largest cluster value; generally represents large adults
w.i <- c(rep(1, nind), rep(0,nz)) #initialize burrows as real/not

### Run the model 
jd.test2 <- list(nind= nind, nz = nz, L = L, Bx = Bx, x=x,y=y,z=z, Nclust = Nclust, q=q, o=o, v=v, v2 = v2, q= q,
                clust= clust, onesRepNclust = rep(1,Nclust))
ji.test2 <-list(w = w.i, betaClust =c(1.6, .8, 1.25, 1.2, 1.3), pClust = c(.4, .25, .2, .1,.05 ),
                sigma.beta = .025, sigma.int = 1, sigma.gamma = -1.5,
                shape = c(52, 9, 65, 55, 60), z.beta = -0.016) 
jp.test2 <- c("D","Dt", "N","Nt", "betaClust", "pClust","shape","sigma.int", "sigma.beta", "p.online", 
              "b.point", "Occ", "sigma.gamma", "juvi1", "juvi2", "juvi3", "z.beta", "o.int")

Foo2 <- autorun.jags(modelstring.Vegetation, data = jd.test2, inits = list(ji.test2, ji.test2, ji.test2), monitor= jp.test2, 
                     startburnin = 10000, startsample = 40000, adapt = 1000, method = 'parallel', n.chain = 3, max.time="8h", 
                     psrf.target = 1.1,silent.jags = FALSE, summarise = TRUE)
Covmod2 <- summary(Foo2)



results <-  list()
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
