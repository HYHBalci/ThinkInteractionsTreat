#' # Think Interactions, Logistic
#' Script: Mark van de Wiel, mark.vdwiel@amsterdamumc.nl
#' Software implements a novel linked shrinkage model for two-way interactions, 
#'  as presented in manuscript:
#' "Linked shrinkage to improve estimation of interaction effects in regression models.
#' In addition, it implements the Shapley values derived in the manuscript.
#' Below's script demonstrates the use of the software for logistic regression. 
#'
#'
#' ## Data preparation
pmt <- proc.time()
setwd("C:\\ExternData\\Helius\\Interactions\\Demo")
show(load("datasynth_Chol.Rdata"))
load("datasynth_Chol.Rdata")

dim(datasynth)
nsam <- nrow(datasynth)

#' Random subset as training
ntrain <- 1000
subset <- sample(1:nsam, ntrain)
data1 <- datasynth[subset,]

#' Creates binary outcome
Y1 <- (sign(as.numeric(data1[,1]))+1)/2
X1 <- as.matrix(data1[,-1])
dim(X1)

#' Source the functions Rstan code and other functions
source('C:/Synchr/Rscripts/SimpleScripts/LinRegr/Interactions/auxiliarycodeRstan.R')
source('C:/Synchr/Rscripts/SimpleScripts/LinRegr/Interactions/auxiliarycodeOther.R')

#' ## New model: Bayint
#' load RStan
library(rstan)

#' Compiles the code for the Bayint model (logistic). Can take a while (e.g. 30 sec).
pmt <- proc.time()
stan_binint_brl_HalfCauchy <- stan_model(model_name = "stan_binint_brl_HalfCauchy", 
                                         model_code = code_binint_brl_HalfCauchy)
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

#' Prepares data and parameter input for Stan.
dat <- list(y = Y1, X = Xmain, Xint = Xint, n = nrow(Xmain), p = p, q = q, g=1:p, g1=g1, g2=g2)

#' Sampling from the model
niter = 10000; nwarmup = 1000 
Bayint <- suppressWarnings(sampling(stan_binint_brl_HalfCauchy, data = dat, iter = niter, 
                                    warmup = nwarmup, thin = 10, chain = 1))

#' Summary of the fit
print(summary(Bayint)$summary)

#' Posteriors of main effect coefficients
betapost <- extract(Bayint, "beta")$beta 
dim(betapost)

#' Posteriors of interaction effect coefficients
betapostint <- extract(Bayint, "betaint")$betaint
dim(betapostint)
