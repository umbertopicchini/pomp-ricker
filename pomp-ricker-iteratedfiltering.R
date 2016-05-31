# Iterated filtering (IF2) estimation of a stochastic Ricker model.
#
# The message is that here IF2 works very well but can fail if starting values are very far from the MLE
# (see results denoted as m3 below), whereas synthetic likelihoods works well also with badly initialized starting values
# (see pomp_ricker-synlik.R).

# This script assumes that the R "pomp" package is already installed.
# For efficiency some construct require compilation.
# If you are on a Windows machine you probably need to install Rtools.
# On Mac you might need Xcode.
# On Linux you should be already good to go.

# Umberto Picchini, 2016.


library("pomp")

state <- Csnippet("
     e = rnorm(0, sigma);
     N = r * N * exp(-N + e);
  ")

rmeas <- Csnippet("
      y = rpois(phi * N);
      ")

dmeas <- Csnippet("
      lik = dpois(y, phi * N, give_log);
      ")

log.trans <- Csnippet("
     Tr = log(r);
     Tsigma = log(sigma);
     Tphi = log(phi);
     TN_0 = log(N_0);
     ")

exp.trans <- Csnippet("
     Tr = exp(r);
     Tsigma = exp(sigma);
     Tphi = exp(phi);
     TN_0 = exp(N_0);
     ")


ricker1 <- pomp(data = data.frame(time = seq(0, 50, by = 1), y = NA),
               rprocess = discrete.time.sim(step.fun = state,
                                            delta.t = 1),
               rmeasure = rmeas,
               dmeasure = dmeas,
               toEstimationScale = log.trans,
               fromEstimationScale = exp.trans,
               paramnames = c("r", "sigma", "phi", "N.0", "e.0"),
               statenames = c("N", "e"), times = "time", t0 = 0,
               params = c(r = exp(3.8), sigma = 0.3, phi = 10, N.0 = 7, e.0 = 0))

ricker1 <- simulate(ricker1, seed = 73691676L)
plot(ricker1)

# redefine the model to work with r=exp(3.799), instead of r=exp(3.8),
ricker2 <- pomp(data = data.frame(time = seq(0, 50, by = 1), y = NA),
                rprocess = discrete.time.sim(step.fun = state,
                                             delta.t = 1),
                rmeasure = rmeas,
                dmeasure = dmeas,
                toEstimationScale = log.trans,
                fromEstimationScale = exp.trans,
                paramnames = c("r", "sigma", "phi", "N.0", "e.0"),
                statenames = c("N", "e"), times = "time", t0 = 0,
                params = c(r = exp(3.799), sigma = 0.3, phi = 10, N.0 = 7, e.0 = 0))

# simulate a new model with the same seed for pseudo-random numbers generation
ricker2 <- simulate(ricker2, seed = 73691676L)
# plot and compare with the previous plot of ricker1, notice the difference in the y path
plot(ricker2)

# consider again ricker1 and perform maximum likelihood estimation via iterated filtering
ricker1 <- simulate(ricker1, seed = 73691676L)
# set parameters starting values (also for constants that need no estimation)
thetastart <- c(r=exp(2.5) , sigma=0.3, phi=10,N.0=7,e.0=0)
# set the seed for pseudo-random numbers (for reproducibility)
set.seed(635363)
estpars <- c("r")  # we only conduct inference for r
theta.guess <- thetastart
# iterated filtering
m1 <- mif2(
             ricker1,
             Nmif=1500,
             start=theta.guess,
             transform=TRUE,
             rw.sd=rw.sd(r=0.1),
             cooling.fraction.50=0.9,
             Np=1000
             )
# restart from the last iteration and perform some more iterations
m1 <- continue(m1,Nmif=200,cooling.fraction=0.7)
m1 <- continue(m1,Nmif=200,cooling.fraction=0.4)
m1 <- continue(m1,Nmif=200,cooling.fraction=0.3)
m1 <- continue(m1,Nmif=200,cooling.fraction=0.2)
plot(m1)  # the algorithm nicely converges close to the true value r=44.7
m1 

# now let's try to see what happens if we consider data generated with log(r)=3.8 and sigma=0.01
ricker3 <- pomp(data = data.frame(time = seq(0, 50, by = 1), y = NA),
                rprocess = discrete.time.sim(step.fun = state,
                                             delta.t = 1),
                rmeasure = rmeas,
                dmeasure = dmeas,
                toEstimationScale = log.trans,
                fromEstimationScale = exp.trans,
                paramnames = c("r", "sigma", "phi", "N.0", "e.0"),
                statenames = c("N", "e"), times = "time", t0 = 0,
                params = c(r = exp(3.8), sigma = 0.01, phi = 10, N.0 = 7, e.0 = 0))

ricker3 <- simulate(ricker3, seed = 73691676L)

# perform iterated filtering on ricker3
# iterated filtering
estpars <- c("r")  # we only conduct inference for r
# set parameters starting values (also for constants that need no estimation)
thetastart3 <- c(r=exp(2.5) , sigma=0.01, phi=10,N.0=7,e.0=0)
theta.guess3 <- thetastart3
set.seed(635363)
m3 <- mif2(
  ricker3,
  Nmif=1500,
  start=theta.guess3,
  transform=TRUE,
  rw.sd=rw.sd(r=0.1),
  cooling.fraction.50=0.9,
  Np=1000
)

m3 <- continue(m3,Nmif=200,cooling.fraction=0.7)
m3 <- continue(m3,Nmif=200,cooling.fraction=0.4)
m3 <- continue(m3,Nmif=200,cooling.fraction=0.3)
m3 <- continue(m3,Nmif=200,cooling.fraction=0.2)
plot(m3) # the algorithm fails in estimating r correctly
m3

# let's try to see what happens if we start from a more favourable value for r (start at r=exp(3.5) instead of r=exp(2.5))
thetastart4 <- c(r=exp(3.5) , sigma=0.01, phi=10,N.0=7,e.0=0)
theta.guess4 <- thetastart4
set.seed(635363)
m4 <- mif2(
  ricker3,
  Nmif=1500,
  start=theta.guess4,
  transform=TRUE,
  rw.sd=rw.sd(r=0.1),
  cooling.fraction.50=0.9,
  Np=1000
)

m4 <- continue(m4,Nmif=200,cooling.fraction=0.7)
m4 <- continue(m4,Nmif=200,cooling.fraction=0.4)
m4 <- continue(m4,Nmif=200,cooling.fraction=0.3)
m4 <- continue(m4,Nmif=200,cooling.fraction=0.2)
plot(m4) # ok now it works ok!
m4


