cholesterol <- TRUE
library(multiridge)
library(mgcv)
library(glmnet)
library(shrinkage)
library(glinternet)
library(rstan)

source('C:/Synchr/Rscripts/SimpleScripts/LinRegr/Interactions/auxiliarycode.R')

if(cholesterol) 
  {
  setwd("C:\\ExternData\\Helius\\Interactions\\Chol\\synth") 
  load("datasynth_Chol.Rdata")
  } else {
  setwd("C:\\ExternData\\Helius\\Interactions\\Sbp\\synth")
    load("datasynth_Sbp.Rdata")
  }
getwd()


#makes 100 data chunks of size 1000
td1000_100 <- makeTraining(data=datasynth,chunksize=1000,nfit=100,nrepeat=5)

#determine columns corresponding to main effects and interactions
cn <- colnames(datasynth[,-1])
mainint <- mapint(cn)
mains <- mainint[[1]]
ncov <- length(cn)
ncovmain <- length(mains)
ncovint <- ncov-ncovmain
c(ncov,ncovmain,ncovint)

#prepare fitting Bayesian models in RStan
pmt <- proc.time()
stan_lin_brl_HalfCauchy <- stan_model(model_name = "stan_lin_brl_HalfCauchy", model_code = code_lin_brl_HalfCauchy)
proc.time()-pmt

pmt <- proc.time()
stan_linint_brl_HalfCauchy <- stan_model(model_name = "stan_linint_brl_HalfCauchy", model_code = code_linint_brl_HalfCauchy)
proc.time()-pmt

#fit OLS, ridge2, Baylcc, Bayint
nset <- 25
trainind <- td1000_100
niter <- 25000; nwarm <- 5000 #mcmcsetting when n<=1000
fitbay <- TRUE
fitbayint <- TRUE
allmodels<- list()
pmttot <- 0
pmt <-proc.time()
for(setind in 1:nset){
  #setind<-2
  print(setind)
  data1 <- datasynth[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  
  #OLS
  lmsub <- lm(resp ~ ., data=data1)
  #summary(lmsub)
  
  #ridge2, using mgcv
  whmain <- 1:ncovmain
  whint <- (1:ncov)[-whmain]
  diag1 <- rep(0,ncov);diag1[whmain] <- 1;diag2 <- rep(0,ncov);diag2[whint] <- 1;
  PP2 <- list(X1=list(diag(diag1),diag(diag2)),sp=-1)
  pred2 <- gam(Y1 ~ 1+X1,family="gaussian",paraPen=PP2,method="ML",scale=-1)
  
  #STAN, Bayloc
  dat <- list(y = Y1, X = X1, n = nrow(X1), p = ncov, g=1:ncov)
  if(fitbay) Bay <- sampling(stan_lin_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarm, thin = 10, chain = 1) else Bay <- "notfit"
  
  #STAN, Bayint (linked shrinkage) 
  Xmain <- X1[,mains]
  Xint <- X1[,-mains]
  g1 <- mainint[[2]][1,]
  g2 <- mainint[[2]][2,]
  
  dat <- list(y = Y1, X = Xmain, Xint = Xint, n = nrow(Xmain), p = ncovmain, q=ncovint, g=1:ncovmain, g1=g1, g2=g2)
  pmtbayint <- proc.time()
  if(fitbayint) Bayint <- sampling(stan_linint_brl_HalfCauchy, data = dat, iter = niter, warmup = nwarm, thin = 10, chain = 1) else Bayint <- "nofit"
  pmtbayint2 <- proc.time() - pmtbayint
  pmttot <- pmttot + pmtbayint2
  print(pmttot)
  print(summary(Bayint)$summary[1:3,1:3])
  models <- list(ols=lmsub,ridge2=pred2,bayloc=Bay,bayint = Bayint)
  allmodels <- c(allmodels,list(models))
  timeelaps <- proc.time()-pmt
  save(timeelaps,allmodels,file="allmodels.Rdata")
}


#lasso variations
catind <- c(3,6)
mat <- c() #interactions for hlasso
for(i in 1:(ncovmain-1)) for(j in (i+1):ncovmain) {
  cond <- (i >= catind[1]) & (i <= catind[2]) & (j >= catind[1]) & (j <= catind[2])
  if(!cond) mat <- rbind(mat,c(i,j))
}


nset <- 25
trainind <- td1000_100
allmodels_lasso<- list()
allmodels_hlasso<- list()

pmt <-proc.time()
for(setind in 1:nset){
  print(setind)
  data1 <- datasynth[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  ncov <- ncol(X1)
  whmain <- 1:ncovmain
  whint <- (1:ncov)[-whmain]
  
  #lasso with unpenalized main effects
  pf <- c(rep(0,length(whmain)),rep(1,length(whint)))
  cvlasso1 <- cv.glmnet(X1,Y1, penalty.factor=pf)
  print(cvlasso1)
  lasso1 <- glmnet(X1,Y1, penalty.factor=pf)
  lasso1$lambda.min <- cvlasso1$lambda.min
  allmodels_lasso <- c(allmodels_lasso,list(lasso1=lasso1))
  
  
  Xmain <- X1[,mains]
  ncovmain <- ncol(Xmain)
  numLevels <- rep(1, ncovmain) 
  #hierarchical lasso
  hlasso <- glinternet.cv(Xmain, Y1, numLevels=rep(1,ncovmain),interactionPairs=mat)
  allmodels_hlasso <- c(allmodels_hlasso,list(hlasso=hlasso))
  timeelaps <- proc.time()-pmt
  save(timeelaps,allmodels_lasso,allmodels_hlasso,file="allmodels_lasso.Rdata")
}

  
#Two-step 
nset <- 25
trainind <- td1000_100
allmodels_2step <- list()
pmt <-proc.time()
for(setind in 1:nset){
  #setind<-25
  print(setind)
  data1 <- datasynth[trainind[[setind]],]
  Y1 <- as.numeric(data1[,1])
  X1 <- as.matrix(data1[,-1])
  cn <- colnames(X1)
  mainint <- mapint(cn)
  mains <- mainint[[1]]
  Xmain <- X1[,mains]
  ncovmain <- ncol(Xmain)
  
  #model with main effects only (step 1)
  lmmain <- lm(Y1 ~ Xmain)
  summ <- summary(lmmain)
  pvalsmain <- summ$coefficients[-1,4]
  whsig <- which(pvalsmain<=0.05)
  datammsel <- data.frame(Xmain[,whsig,drop=FALSE])
  
  #include interactions of significant main effects only
  mmsel <- model.matrix(~(.)^2, data=datammsel)[,-1,drop=FALSE]
  matchcn <- match(colnames(mmsel), cn)
  whna <- which(is.na(matchcn))
  if(length(whna)>0) matchcn <- matchcn[-whna]
  
  #step 2: fit model with main effects and selected interactions
  lmtwostep <- lm(Y1 ~ . , data = data.frame(X1[,matchcn,drop=FALSE]))
  lmtwostep$matchcn <- matchcn
  allmodels_2step <- c(allmodels_2step,list(lmtwostep))
  timeelaps <- proc.time()-pmt
  save(timeelaps,allmodels_2step,file="allmodels_2step.Rdata")
}



load("allmodels.Rdata")
load("allmodels_lasso.Rdata")
load("allmodels_2step.Rdata")


#Retrieve coefficients from all models and subsets

#function to retrieve coefs from hlasso
coefshlasso <- function(cv,cndat,mat){
  #cv=allmodels_hlasso[[i]];cndat=colnames(helius4int[,-1])
  betahat <- rep(0, length(cndat)+1)   
  names(betahat) <- c("(Intercept)", cndat)
  betahat[1] <- cv$betahat[[1]][1]
  betahat[1 + coef(cv)$mainEffects$cont] <- unlist(coef(cv)$mainEffectsCoef$cont)
  interactions <- coef(cv)$interactions$contcont
  wh <- match(asplit(interactions,1),asplit(mat,1))
  betahat[1 + ncovmain + wh] <- unlist(coef(cv)$interactionsCoef$contcont)
  return(betahat)
}


#true (OLS coef from master set)
lmmaster <- lm(resp ~ ., data = datasynth)
true <- coef(lmmaster)[-1]

#retrieve coefficients
ncoef <- ncol(datasynth)-1
cn <- colnames(datasynth)[-1]
allcoef <- list()
for(i in 1:length(allmodels)){
  #i <- 1
  models <- allmodels[[i]]
  coefols <- coef(models$ols)[-1] #exclude intercept
  coefridge2 <- coef(models$ridge2)[-1] #exclude intercept
  coefbay <- summary(models$bayloc)$summary[3:(2+ncoef),1]  #exclude intercept and hyperparameter 
  coefbayint <- summary(models$bayint)$summary[4:(3+ncoef),1] #exclude intercept and 2 hyperparameters
   lasso1 <- allmodels_lasso[[i]]
  coeflas <- as.numeric(coef(lasso1,s=lasso1$lambda.min))[-1]
  coefhlas <- coefshlasso(cv=allmodels_hlasso[[i]],cndat=cn,mat=mat)[-1]
  mod2step <- allmodels_2step[[i]]
  whsel <- mod2step$matchcn
  coef2step <- rep(0,ncoef)
  coef2step[whsel] <- coef(mod2step)[-1]
  allcoefi <- cbind(true, coefols,coefridge2,coefbay,coefbayint,coeflas,coefhlas,coef2step)
  allcoefi
  allcoef <- c(allcoef,list(allcoefi))
}



#compute root MSEs
se <- lapply(allcoef, function(aci) {
  sqer <- (aci[,-1] - aci[,1])^2
  colnames(sqer) <- colnames(allcoef[[1]])[-1]
  return(sqer)
})

rmsek <- function(k){
  #k<-1
  msek <- sapply(se,function(mat) mat[,k] )
  rmsek <- sqrt(apply(msek,1,mean))
  return(rmsek)
}
rmses <- sapply(1:ncol(se[[1]]),rmsek)
colnames(rmses) <- colnames(allcoef[[1]])[-1]

#retrieve significant effects from master set (for plotting)
sumwhole <- summary(lmmaster)
pvals <- sumwhole$coefficients[-1,4]
pvordm <- order(pvals[mains])
pvordi <- order(pvals[-mains])
nsigm <- sum(pvals[mains] <= 0.01)
nsigm
nsigi <- sum(pvals[-mains] <= 0.01)
nsigi


#plotting
#four plots IN MANUSCRIPT
if(cholesterol) pdf(file="rMSE_interact_chol_synth.pdf",width=14,height=10) else pdf(file="rMSE_interact_sbp_synth.pdf",width=14,height=10)
par(mfrow=c(2,2))
toplot1 <- c(1:3,4)
toplot2 <- c(5:7,4)
lwds <- c(rep(2,3),3)
cols1 <- c(2:4,1)
cols2 <- c(5:7,1)
for(comp in c(1,2)){
  if(comp==1){toplot<-toplot1;cols<-cols1} else {toplot<-toplot2;cols<-cols2}
  rmsesord <- rbind(rmses[mains[pvordm],],rmses[-mains,][pvordi,])
  allrMSEsel <- rmsesord[,toplot]
  rncov <- rownames(allrMSEsel)
  ncov <- nrow(allrMSEsel)
  maxrmse <- max(allrMSEsel)
  maxrmse
  par(las=2,xaxt="n",mar=c(8,4,2,1))
  if(comp==1) mn = "rMSEs" else mn = ""
  plot(1:ncov,allrMSEsel[,1],col=cols[1],type="o",lwd=2,xlab="",ylab="rMSE",xlim = c(1,ncov),
       ylim=c(0,0.17),main=mn, cex=0.8)
  for(i in 1:ncol(allrMSEsel)) points(1:ncov,allrMSEsel[,i],col=cols[i],pch=i,type="o",lwd=lwds[i],
                                      cex=0.8)
   par(las=2,xaxt="s")
  if(comp==1) plotlegend <- c("OLS","ridge2","Bayloc","Bayint") else  plotlegend <- c("lassoint","hlasso","2step","Bayint")
  legend(0.88*ncov,0.17,legend=plotlegend,col=cols,
         lwd=lwds,pch=toplot,cex=0.8)
  abline(v=mains[length(mains)] + 0.5,lwd=3)
  abline(v= nsigm + 0.5,lwd=1)
  abline(v= length(mains) + nsigi + 0.5,lwd=1)
  
  
  
  #Zooming in
  rmsesord <- rbind(rmses[mains[pvordm],],rmses[-mains,][pvordi,])[(1:(ncovmain + 2*nsigi)),]
  allrMSEsel <- rmsesord[,toplot]
  rncov <- rownames(allrMSEsel)
  ncov <- nrow(allrMSEsel)
  maxrmse <- max(allrMSEsel)
  maxrmse
  par(las=2,xaxt="n",mar=c(8,4,2,1))
  if(comp==1) mn = "rMSEs, zoom" else mn = ""
  plot(1:ncov,allrMSEsel[,1],col=cols[1],type="o",lwd=2,xlab="",ylab="rMSE",xlim = c(1,ncov),
       ylim=c(0,0.17),main=mn, cex=0.8)
  for(i in 1:ncol(allrMSEsel)) points(1:ncov,allrMSEsel[,i],col=cols[i],pch=i,type="o",lwd=lwds[i],
                                      cex=0.8)
  par(las=2,xaxt="s")
  axis(1,at=1:ncov,rncov)

  if(comp==1) plotlegend <- c("OLS","ridge2","Bayloc","Bayint") else  plotlegend <- c("lassoint","hlasso","2step","Bayint")
  legend(0.88*ncov,0.17,legend=plotlegend,col=cols,lwd=lwds,
         pch=1:7,cex=0.8)
  abline(v=mains[length(mains)] + 0.5,lwd=3)
  abline(v= nsigm + 0.5,lwd=1)
  abline(v= length(mains) + nsigi + 0.5,lwd=1)
}
dev.off()





