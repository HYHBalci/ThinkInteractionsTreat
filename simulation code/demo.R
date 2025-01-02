source('auxiliarycodeRstan.R')
source('auxiliarycodeOther.R')
source('simulation code/simulations.R')
p_main <- 5
p_noise_main <- 4
p <- p_main + p_noise_main
data <- simulate_data(n_samples = 100, p_main = p_main, p_noise_main = p_noise_main, interaction = TRUE, treatment = TRUE,seed = 123)

library(rstan)
data$main_effects

data1 <- data$simulated_data
Y1 <- as.numeric(data1[,1])
A <- as.numeric(data1[,47])
X1 <- as.matrix(data1[,-c(1, ncol(data1))])
dim(X1)
#' load RStan
library(rstan)

#' Compiles the code for the Bayint model. Can take a while (e.g. 30 sec).
pmt <- proc.time()
stan_linint_brl_HalfCauchy_treat <- stan_model(model_name = "stan_linint_brl_HalfCauchy_treat", 
                                         model_code = code_linint_brl_HalfCauchy_treat)
proc.time()-pmt

#' Create submatrices with main effects and interactions
cn <- colnames(X1)
mainint <- mapint(cn) 
mains <- mainint[[1]]
Xmain <- X1[,mains]
Xint <- X1[,-mains]

#' Nr of main effects
p <- ncol(Xmain)
p

#Nr of interactions
q <- ncol(Xint)
q

#' Indexing for linking the interactions with main effects. 
g1 <- mainint[[2]][1,] #contains main effect indices for 1rst term in interaction
g2 <- mainint[[2]][2,] #contains main effect indices for 2nd term in interaction
g1

#' Prepares data and parameter input for Stan.
dat <- list(y = Y1, X = Xmain, Xint = Xint, a = A, n = nrow(Xmain), p = p, q = q, g=1:p, g1=g1, g2=g2)

#' Sampling from the model
niter = 25000; nwarmup = 5000 

#' NOTE 1: RECOMMENDED TO USE 25000, 5000 FOR FINAL RESULTS 
Bayint <- suppressWarnings(sampling(stan_linint_brl_HalfCauchy_treat, data = dat, iter = niter, warmup = nwarmup, thin = 10, chain = 1))

#' Summary of the fit
print(summary(Bayint)$summary)

# Extract posterior samples
posterior_samples <- as.data.frame(Bayint)

# Example: Plot posterior distribution for "alpha"
hist(posterior_samples$alpha, breaks = 30, main = "Posterior of alpha", xlab = "alpha")

# Example: Plot posterior distribution for "gamma"
hist(posterior_samples$gamma, breaks = 30, main = "Posterior of gamma", xlab = "gamma")