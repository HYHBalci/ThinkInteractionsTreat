#' # Think Interactions
#' Script: Mark van de Wiel, mark.vdwiel@amsterdamumc.nl
#' Software implements a novel linked shrinkage model for two-way interactions, 
#'  as presented in manuscript:
#' "Linked shrinkage to improve estimation of interaction effects in regression models.
#' In addition, it implements the Shapley values derived in the manuscript.
#' Below's script demonstrates the use of the software. 
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
Y1 <- as.numeric(data1[,1])
X1 <- as.matrix(data1[,-1])
dim(X1)

#' Source the functions Rstan code and other functions
source('C:/Synchr/Rscripts/SimpleScripts/LinRegr/Interactions/auxiliarycodeRstan.R')
source('C:/Synchr/Rscripts/SimpleScripts/LinRegr/Interactions/auxiliarycodeOther.R')

#' ## New model: Bayint
#' load RStan
library(rstan)

#' Compiles the code for the Bayint model. Can take a while (e.g. 30 sec).
pmt <- proc.time()
stan_linint_brl_HalfCauchy <- stan_model(model_name = "stan_linint_brl_HalfCauchy", 
                                         model_code = code_linint_brl_HalfCauchy)
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
dat <- list(y = Y1, X = Xmain, Xint = Xint, n = nrow(Xmain), p = p, q = q, g=1:p, g1=g1, g2=g2)

#' Sampling from the model
niter = 10000; nwarmup = 1000 

#' NOTE 1: RECOMMENDED TO USE 25000, 5000 FOR FINAL RESULTS 
Bayint <- suppressWarnings(sampling(stan_linint_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarmup, thin = 10, chain = 1))

#' Summary of the fit
print(summary(Bayint)$summary)

#' Posteriors of main effect coefficients
betapost <- extract(Bayint, "beta")$beta 
dim(betapost)

#' Posteriors of interaction effect coefficients
betapostint <- extract(Bayint, "betaint")$betaint
dim(betapostint)

#' ## Variable importance 
ntest <- 100
subsettest <- sample(1:nsam, ntest)
datatest <- datasynth[subsettest,]
Ytest <- as.numeric(datatest[,1])
Xtest <- as.matrix(datatest[,-1])

#' Personalized unit change posterior of "age" for all test samples
postallage <- postallf("age", Xtest, betapostb=betapost, betaposti=betapostint)

#' Plot posterior intervals of personalized unit change of age for all test samples
gender <- datatest$genderf
ord <- order(datatest$age)
genderord <- gender[ord]
paord <- postallage[,ord]
par(mar=c(1,2,2,1))
plot(paord[2,],pch=16+genderord, col=genderord+3, ylim =c(-0.4,.7),xlab="",ylab="")
segments(1:ntest,paord[1,],1:ntest,paord[3,], col=genderord+3)
legend(85,-0.25,legend=c("male", "female"),pch=c(15,17),col=c(2,4))

#' ### Shapley values 
#' 
#' Needs means and (co-)variances, as computed from 
EXX <- cov(X1)
EX <- apply(X1,2,mean)

#'  Covariate for which Shapley values are computed. May be a name (should 
#'  match with that in the data set) or an index
cov <- "age"

#' Compute Shapley values and intervals. Three components: 1: Total; 2: 
#' Contribution main effect; 3: contribution interactions. 
shapcov <- shapley(cov, Xtest ,EXX, EX,betapost,betaposti) 

#' Plot Shapleys and 95% cred int for test samples ordered w.r.t. covariate values
par(mar=c(2,2,2,1),mfrow=c(1,3))
for(k in 1:3){
shapcovtot <- shapcov[[k]] #first component is total shapley value
whcov <- which(colnames(Xtest)==cov)
ord <- order(Xtest[,whcov])
covvalue <- Xtest[,whcov]
shord <- shapcovtot[ord,]
par(mar=c(2,2,2,1))
if(k==1) lab <- "Shapley values, Total"
if(k==2) lab <- "Shapley values, Main"
if(k==3) lab <- "Shapley values, Inter"
plot(shord[,2],pch=15, col=1,xlab="",ylim=c(-0.7,0.7),ylab="",main=lab)
segments(1:ntest,shord[,1],1:ntest,shord[,3], col=1)
abline(h=0)
}

#' Compute global variable importances from Shapleys across large test sample (I_j)
ntestl <- 1000
subsettestl <- sample(1:nsam, ntestl)
datatestl <- datasynth[subsettestl,]
Xtestl <- as.matrix(datatest[,-1])

#' Plot global Shaps for all covariates
allshap <- shapleys(X1, Xtestl, betapostb,betaposti,grouped = list(etn=3:6))
allshapsum <- sapply(1:length(allshap), function(i) return(allshap[[i]]$tot[,2]))
colnames(allshapsum) <-  names(allshap)
par(mar=c(3,2,2,1))
boxplot(allshapsum, main="Shapleys Total",ylim = c(-0.7,0.7))


#' ## Variations to Bayint
#' Bay0int: no shrinkage of main effects (vague prior)
stan_lin0int_brl_HalfCauchy <- stan_model(model_name = "stan_lin0int_brl_HalfCauchy", model_code = code_lin0int_brl_HalfCauchy)
Bay0int <- suppressWarnings(sampling(stan_lin0int_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarmup, thin = 10, chain = 1))

#' Bayintadd: additive model for shrinkage of interaction parameters
stan_linintadd_brl_HalfCauchy <- stan_model(model_name = "stan_linintadd_brl_HalfCauchy", model_code = code_linintadd_brl_HalfCauchy)
Bayintadd <- suppressWarnings(sampling(stan_linintadd_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarmup,  thin = 10, chain = 1))

#' Bayint*: no adaptive global shrinkage parameter (tau2int=1) for interaction parameters
stan_linint2_brl_HalfCauchy <- stan_model(model_name = "stan_linint2_brl_HalfCauchy", model_code = code_linint2_brl_HalfCauchy)
dat <- list(y = Y1, X = Xmain, Xint = Xint, n = nrow(Xmain), p = p, q=q, tau2int=1, g=1:p, g1=g1, g2=g2)
Bayintstar <- suppressWarnings(sampling(stan_linint2_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarmup, thin = 10, chain = 1))

#' ## Alternative methods
#' ### OLS
lmsub <- lm(resp ~ ., data=data1)
summary(lmsub)

#' Bayloc: Standard local shrinkage model
stan_lin_brl_HalfCauchy <- stan_model(model_name = "stan_lin_brl_HalfCauchy", model_code = code_lin_brl_HalfCauchy)
nc <- ncol(X1)
dat <- list(y = Y1, X = X1, n = nrow(X1), p = nc, g=1:nc)
Bayloc <- suppressWarnings(sampling(stan_lin_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarmup, thin = 10, chain = 1))

#' ### ridge2
#' load mgcv, which is used to estimate 2 ridge penalties automatically
library(mgcv)
ncov <- ncol(X1)

#' Retrieve number of main effects from column names X1
cn <- colnames(X1)
mainint <- mapint(cn) 
mains <- mainint[[1]]
ncovmain <- length(mains)

#' Create two penalty matrices, one for main effects, one for interactions
whmain <- 1:ncovmain
whint <- (1:ncov)[-whmain]
diag1 <- rep(0,ncov);diag1[whmain] <- 1;diag2 <- rep(0,ncov);diag2[whint] <- 1;
PP2 <- list(X1=list(diag(diag1),diag(diag2)),sp=-1)

#' Fit with mgcv's gam function
ridge2 <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2,method="ML",scale=-1)
ridge2$sp #estimated penalties
summary(ridge2) #fit summary

#' ### Two-step (add only two-way interactions of significant main effects) 
#' Fit lm to model with main effects only
Xmain <- X1[,mains]
cn <- colnames(X1)
ncovmain <- ncol(Xmain)
lmmain <- lm(Y1 ~ Xmain)

#' Extract significant main effects
summ <- summary(lmmain)
pvalsmain <- summ$coefficients[-1,4]
whsig <- which(pvalsmain<=0.05)
datammsel <- data.frame(Xmain[,whsig,drop=FALSE])

#' Create model matrix with all two-way interactions of significant main effects
mmsel <- model.matrix(~(.)^2, data=datammsel)[,-1,drop=FALSE]
matchcn <- match(colnames(mmsel), cn)
whna <- which(is.na(matchcn))
if(length(whna)>0) matchcn <- matchcn[-whna]

#' Fit final model
lmtwostep <- lm(Y1 ~ . , data = data.frame(X1[,matchcn,drop=FALSE]))

#' Store match of selected columns with those of the original data frame, X1
lmtwostep$matchcn <- matchcn

#' ### lassoint  (lasso penalty on interactions)
library(glmnet)
whmain <- 1:ncovmain
whint <- (1:ncov)[-whmain]
#' penalty factors: 0 for main effects, 1 for interaction terms
pf <- c(rep(0,length(whmain)),rep(1,length(whint)))
cvlasso1 <- cv.glmnet(X1,Y1, penalty.factor=pf)
print(cvlasso1)

#' Fit lasso and store optimal lambda
lasso1 <- glmnet(X1,Y1, penalty.factor=pf)
lasso1$lambda.min <- cvlasso1$lambda.min

#' ### hierarchical lasso 
library(glinternet)
Xmain <- X1[,mains]
ncovmain <- ncol(Xmain)
numLevels <- rep(1, ncovmain)

#' hlasso needs the interactionPairs. As we wish to exclude interactions between 
#' levels of the categorical variable, which represent dummies 3 to 6, this requires
#' some programming. Below's code creates all interactionPairs between covariates 1, .., 14
#' except those between 3, ..., 6.
#' 
catind <- c(3,6)
mat <- c() #interactions for hlasso
for(i in 1:(ncovmain-1)) for(j in (i+1):ncovmain) {
  cond <- (i >= catind[1]) & (i <= catind[2]) & (j >= catind[1]) & (j <= catind[2])
  if(!cond) mat <- rbind(mat,c(i,j))
}
print(mat)

#' Fit and CV hlasso
hlassocv <- glinternet.cv(Xmain, Y1, numLevels=rep(1,ncovmain),interactionPairs=mat)
proc.time()-pmt

#' ## Retrieving coefficients
#' function to retrieve coefs from hlasso; a bit more tedious than for the 
#' methods
coefshlasso <- function(cv,cndat,mat){
  betahat <- rep(0, length(cndat)+1)   
  names(betahat) <- c("(Intercept)", cndat)
  betahat[1] <- cv$betahat[[1]][1]
  betahat[1 + coef(cv)$mainEffects$cont] <- unlist(coef(cv)$mainEffectsCoef$cont)
  interactions <- coef(cv)$interactions$contcont
  wh <- match(asplit(interactions,1),asplit(mat,1))
  betahat[1 + ncovmain + wh] <- unlist(coef(cv)$interactionsCoef$contcont)
  return(betahat)
}

#' Estimated coefficients for all methods (posterior means for Bayesian estimates)
ncoef <- ncol(X1)
coefBayint <- summary(Bayint)$summary[4:(3+ncoef),1] #exclude intercept and 2 hyperparameters
coefBay0int <- summary(Bay0int)$summary[3:(2+ncoef),1];
coefBayintadd <- summary(Bayintadd)$summary[4:(3+ncoef),1]
coefBayintstar <- summary(Bayintstar)$summary[3:(2+ncoef),1]
coefols <- coef(lmsub)[-1] #exclude intercept
coefridge2 <- coef(ridge2)[-1] #exclude intercept
coefBayloc <- summary(Bayloc)$summary[3:(2+ncoef),1] #exclude intercept and hyperparameter 
coeflassoint <- as.numeric(coef(lasso1,s=lasso1$lambda.min))[-1]
coefhlas <- coefshlasso(cv=hlassocv,cndat=colnames(X1),mat=mat)[-1]
whsel <- lmtwostep$matchcn
coef2step <- rep(0,ncoef)
coef2step[whsel] <- coef(lmtwostep)[-1]
allcoef <- data.frame(coefBayint,coefBay0int,coefBayintadd,coefBayintstar,coefols,
                      coefridge2, coefBayloc,coeflassoint,coefhlas,coef2step)
print(allcoef)


#' ## Logistic regression



