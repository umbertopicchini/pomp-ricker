# Particle marginal estimation of a stochastic Ricker model.
# Umberto Picchini, 2016
# www.maths.lth.se/matstat/staff/umberto/

# We show that, unless the particle MCMC algorithm is correctly initialised with a good starting value
# for the parameter r, the algorithm fails. Compare with pomp_ricker-synlik.R where the synthetic likelihoods
# algorithm is able to approach true value for r even when the starting value is far from the  true one.

# This script assumes that the R "pomp" package is already installed.
# For efficiency some construct require compilation.
# If you are on a Windows machine you probably need to install Rtools.
# On Mac you might need Xcode.
# On Linux you should be already good to go.

# load pomp
library("pomp")


# define states equation and observations model for a stochastic Ricker model
ricker.sim  <- "e = rnorm(0, sigma);N = r * N * exp(-N + e);"  # states equation
ricker.rmeas <- "y = rpois(phi * N);"                          # mechanism to draw from the observations model
ricker.dmeas <- "lik = dpois(y, phi * N, give_log);"           # probability mass function for the observations


# we set a uniform prior on parameter r (below we conduct inference on r only)
densprior <- Csnippet("
    lik = dunif(r,exp(2),exp(5),1);
                      lik = (give_log) ? lik : exp(lik);
                      ")

# useful parameter transformations to ease MCMC sampling
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
       dprior=densprior,toEstimationScale = 
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
thetastart <- c(r=exp(2.5) , sigma=0.3, phi=10,N.0=7,e.0=0)
# set the seed for pseudo-random numbers (for reproducibility)
set.seed(635363)
# use pmcmc to estimate r (1000 particles and 2000 MCMC iterations)
pmcmc1<-pmcmc(myobservedricker, start = thetastart,
               Nmcmc = 2000, Np = 1000, max.fail = Inf,
              proposal=mvn.diag.rw(c(r = 2)))

# plot essential results
plot(pmcmc1)
#extract chains for examination via CODA
library(coda)
library(magrittr)
pmcmc1 %>% conv.rec() -> traces
traceplot(traces[,4])  # traceplot for r
burnin=500
mean(traces[burnin:2000,4])  # compute the posterior mean for r, after some burnin
quantile(traces[burnin:2000,4],probs=c(0.025,0.975))  # 95% posterior bounds for r

# now let's try with data generated with a smaller sigma.
# We expect to obtain bad results as the likelihood approximated via bootstrap filter is supposed to be very poor for small sigma.
# Redefine a "ricker object" using a different sigma. 
ricker2 <- pomp(data = data.frame(time = seq(1, 50, by = 1), y = NA),
          rprocess = discrete.time.sim(step.fun =    
          Csnippet(ricker.sim),delta.t = 1), rmeasure = 
          Csnippet(ricker.rmeas),dmeasure  
          =Csnippet(ricker.dmeas),
          dprior=densprior,toEstimationScale = 
          Csnippet(log.trans),fromEstimationScale = 
          Csnippet(exp.trans),paramnames = c("r", "sigma", 
          "phi", "N.0", "e.0"),statenames = c("N", "e"), times = 
          "time", t0 = 0,params = c(r = exp(3.8), sigma = 0.01, 
           phi = 10,N.0 = 7, e.0 = 0)) 

# generate new observations	(same random seed)
newobservedricker <- simulate(ricker2,seed=73691676L)
# plot observations and latent variables
plot(newobservedricker)
# set parameters starting values (also for constants that need no estimation)
thetastart <- c(r=exp(2.5) , sigma=0.01, phi=10,N.0=7,e.0=0)
# set the seed for pseudo-random numbers (for reproducibility)
set.seed(635363)
# use pmcmc to estimate r 
pmcmc2<-pmcmc(newobservedricker, start = thetastart,
              Nmcmc = 2000, Np = 1000, max.fail = Inf,
              proposal=mvn.diag.rw(c(r = 2)))
# plot essential results
plot(pmcmc2)  # the algorithm fails in sampling the posterior of r

# Finally let's see what happens if we start at a better value for r
thetastart2 <- c(r=exp(3.5) , sigma=0.01, phi=10,N.0=7,e.0=0)
# set the seed for pseudo-random numbers (for reproducibility)
set.seed(635363)
# use pmcmc to estimate r 
pmcmc3<-pmcmc(newobservedricker, start = thetastart2,
              Nmcmc = 2000, Np = 1000, max.fail = Inf,
              proposal=mvn.diag.rw(c(r = 2)))
plot(pmcmc3)  # ok if we start at a value close to the true one the chain converges
