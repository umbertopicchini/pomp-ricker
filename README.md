# pomp-ricker
Examples of how to fit a stochastic Ricker model using the R pomp package.

At the moment provided scripts are:
- pomp_ricker-pmcmc producing exact Bayesian inference using particle MCMC (the particle marginal methodology of Andrieu and Roberts 2009).
- pomp_ricker-synlik producing approximate maximum likelihood estimation using Wood's 2010 synthetic likelihoods.
- pomp_ricker-iteratedfiltering producing approximate maximum likelihood estimation using the improved iterated filtering method of  Ionides et al. 2015

The main purpose of these codes is to show that, unless the parameters starting values are correctly initialized, particle MCMC and iterated filtering may fail when the system is nearly deterministic. For the considered example the synthetic likelihoods methodology is able to converge towards the true value of the parameters even when these are badly initialized.
The code pomp_ricker-synlik also loads the "synlik" package to verify distributional assumptions on the summary statistics.

Our programs have been inspired by codes (or codes snippets) available on the very much recommended resource https://kingaa.github.io/pomp/docs.html, the research article available at https://arxiv.org/abs/1509.00503 and https://cran.r-project.org/web/packages/synlik/vignettes/synlik.html

# Installation
If you are on a Windows machine you probably need to install Rtools.
- On Mac you might need Xcode.
- On Linux you should be already good to go.

References:
- Andrieu, C., & Roberts, G. O. (2009). The pseudo-marginal approach for efficient Monte Carlo computations. The Annals of Statistics, 697-725.
- Wood, S. N. (2010). Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310), 1102-1104.
- Ionides, E. L., Nguyen, D., Atchadé, Y., Stoev, S., & King, A. A. (2015). Inference for dynamic and latent variable models via iterated, perturbed Bayes maps. Proceedings of the National Academy of Sciences, 112(3), 719-724.
