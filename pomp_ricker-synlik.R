# Synthetic likelihoods estimation of a stochastic Ricker model.
# Umberto Picchini, 2016

# This script assumes that the R "pomp" and "synlik" packages are already installed.
# For efficiency some construct require compilation.
# If you are on a Windows machine you probably need to install Rtools.
# On Mac you might need Xcode.
# On Linux you should be already good to go.

# load pomp
library("pomp")


# define states equation and observations model for a stochastic Ricker model
ricker.sim  <- "e = rnorm(0, sigma);N = r * N * exp(-N + e);"  # states equation
ricker.rmeas <- "y = rpois(phi * N);"                          # mechanism to draw from the observation model
ricker.dmeas <- "lik = dpois(y, phi * N, give_log);"           # probability mass function for the observations

# useful parameter transformations to ease omptimization
log.trans <- "Tr = log(r); Tsigma = log(sigma); 
              Tphi = log(phi); TN_0 = log(N_0);"
exp.trans <- "Tr = exp(r);Tsigma = exp(sigma);
              Tphi = exp(phi);TN_0 = exp(N_0);"


# define a "ricker object"             
ricker <- pomp(data = data.frame(time = seq(1, 50, by = 1), y   
       = NA),rprocess = discrete.time.sim(step.fun =    
       Csnippet(ricker.sim),delta.t = 1), rmeasure = 
       Csnippet(ricker.rmeas),dmeasure  
       =Csnippet(ricker.dmeas),
       toEstimationScale = 
       Csnippet(log.trans),fromEstimationScale = 
       Csnippet(exp.trans),paramnames = c("r", "sigma", 
       "phi", "N.0", "e.0"),statenames = c("N", "e"), times = 
       "time", t0 = 0,params = c(r = exp(3.8), sigma = 0.3, 
       phi = 10,N.0 = 7, e.0 = 0)) 

# generate observations	   
myobservedricker <- simulate(ricker,seed=73691676L)
# plot observations and latent variables
plot(myobservedricker)
# set parameters starting values (also for constants that need no estimation)
thetastart <- c(r=exp(2.5) , sigma=1, phi=20,N.0=7,e.0=0)
# set the seed for pseudo-random numbers (for reproducibility)
set.seed(635363)
# the following "probes" are essentially summary statistics that can be applied either to observed and simulated data
plist <- list(probe.marginal("y", ref = obs(myobservedricker), transform = sqrt),
                              probe.acf("y", lags = c(0, 1, 2, 3, 4), transform = sqrt),
                              probe.nlar("y", lags = c(1, 1, 1, 2), powers = c(1, 2, 3, 1),
                              transform = sqrt))

# compute summaries for the observed data
summary.truth <- probe(myobservedricker, probes = plist, nsim = 5000, seed = 1066L)
# compute summaries at parameters starting values
summary.guess <- probe(myobservedricker, params = thetastart, probes = plist, nsim = 5000, seed = 1066L)
# optimize the likelihood approximated via Wood's synthetic likelihoods method
synlikresult <- probe.match(summary.guess, est = c("r", "sigma", "phi"), transform = TRUE,
                      method = "Nelder-Mead", maxit = 10000, seed = 1066L, reltol = 1e-08)
summary(synlikresult)


# now let's try with data generated with a smaller sigma (with this setting particle MCMC has failed)
ricker2 <- pomp(data = data.frame(time = seq(1, 50, by = 1), y = NA),
          rprocess = discrete.time.sim(step.fun =    
          Csnippet(ricker.sim),delta.t = 1), rmeasure = 
          Csnippet(ricker.rmeas),dmeasure  
          =Csnippet(ricker.dmeas),
          toEstimationScale = 
          Csnippet(log.trans),fromEstimationScale = 
          Csnippet(exp.trans),paramnames = c("r", "sigma", 
          "phi", "N.0", "e.0"),statenames = c("N", "e"), times = 
          "time", t0 = 0,params = c(r = exp(3.8), sigma = 0.01, 
           phi = 10,N.0 = 7, e.0 = 0)) 

# generate new observations	(same random seed)
newobservedricker <- simulate(ricker2,seed=73691676L)
plist2 <- list(probe.marginal("y", ref = obs(newobservedricker), transform = sqrt),
              probe.acf("y", lags = c(0, 1, 2, 3, 4), transform = sqrt),
              probe.nlar("y", lags = c(1, 1, 1, 2), powers = c(1, 2, 3, 1),
                         transform = sqrt))
# compute summaries for the observed data
summary.truth2 <- probe(newobservedricker, probes = plist2, nsim = 5000, seed = 1066L)
# compute summaries at parameters starting values
summary.guess2 <- probe(newobservedricker, params = thetastart, probes = plist2, nsim = 5000, seed = 1066L)
# optimize the likelihood approximated via Wood's synthetic likelihoods methos
synlikresult2 <- probe.match(summary.guess2, est = c("r", "sigma", "phi"), transform = TRUE,
                            method = "Nelder-Mead", maxit = 10000, seed = 1066L, reltol = 1e-08)
summary(synlikresult2)

#load the synlik packages to verify gaussianity of summary statistics at some parameter values
library(synlik)
ricker_sl <- synlik(simulator = rickerSimul,
                    summaries = rickerStats,
                    param = c( logR = 3.8, logSigma = log(0.3), logPhi = log(10) ),
                    extraArgs = list("nObs" = 50, "nBurn" = 50)
)
#Now we are ready to simulate a dataset from the object:
ricker_sl@data <- simulate(ricker_sl, nsim = 1, seed = 73691676L)
# Store the data in the object created above 
ricker_sl@extraArgs$obsData <- ricker_sl@data 
# check the approximate normality of statistics or the chi-squared distribution of 
# (s-mu)*inv_sigma*(s-mu)' which is chi-squared(d) (with degrees of freedom d=dim(sobs)).
# See the "Checking the normality assumption and goodness of fit" in the appendix of Wood 2010 for details.
checkNorm(ricker_sl,nsim=1000)



