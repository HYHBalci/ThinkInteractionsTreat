#------------------------------------------------------------------------------#
# Bayesian regression with global and local shrinkage priors
#------------------------------------------------------------------------------#
# global shrinkage with half-Cauchy prior.  Sampling from the half-Cauchy as in 
# Makalic and Schmidt (2016), IEEE Signal Processing Letters:
# tau ~ C+(0,sc) if tau^2 | gam^2 ~ IG(1/2,1/gam^2) and gam^2 ~ IG(1/2,1/sc^2)


##### LINEAR MODEL ###########
# Bayint
code_linint_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> q;
    int<lower=0> n;
    matrix[n,p] X;
    matrix[n,q] Xint;
    real y[n];
    int g[p];
    int g1[q];
    int g2[q];
  }
  parameters {
    real alpha;
    real<lower=0> sigma2;
    real<lower=0.01, upper=1> tau2int;
    vector[p] beta;
    vector[q] betaint;
    real<lower=0>  gam2[p];
    real<lower=0>  tau2[p];
  }
  model {
    sigma2 ~ inv_gamma(1,0.001);
    tau2int ~ uniform(0.01,1);
    alpha ~ normal(0, sqrt(100));
    for (k in 1:p) {
      gam2[k] ~ inv_gamma(0.5, 1);
      tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(sigma2*tau2[g[i]]));
    }
      for (i in 1:q) {
      betaint[i] ~ normal(0, sqrt(sigma2*sqrt(tau2[g1[i]]*tau2[g2[i]])*tau2int));
    }
    y ~ normal(alpha + X*beta + Xint*betaint,sqrt(sigma2));
  }
"

# Bay0int
code_lin0int_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> q;
    int<lower=0> n;
    matrix[n,p] X;
    matrix[n,q] Xint;
    real y[n];
    int g[p];
    int g1[q];
    int g2[q];
  }
  parameters {
    real alpha;
    real<lower=0> sigma2;
    vector[p] beta;
    vector[q] betaint;
    real<lower=0>  gam2[p];
    real<lower=0>  tau2[p];
  }
  model {
    sigma2 ~ inv_gamma(1,0.001);
    alpha ~ normal(0, sqrt(100));
    for (k in 1:p) {
       gam2[k] ~ inv_gamma(0.5, 1);
      tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(100));
    }
      for (i in 1:q) {
      betaint[i] ~ normal(0, sqrt(sigma2*sqrt(tau2[g1[i]]*tau2[g2[i]])));
    }
    y ~ normal(alpha + X*beta + Xint*betaint,sqrt(sigma2));
  }
"

# Bayintadd
code_linintadd_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> q;
    int<lower=0> n;
    matrix[n,p] X;
    matrix[n,q] Xint;
    real y[n];
    int g[p];
    int g1[q];
    int g2[q];
  }
  parameters {
    real alpha;
    real<lower=0> sigma2;
    real<lower=0.01, upper=1> tau2int;
    vector[p] beta;
    vector[q] betaint;
    real<lower=0>  gam2[p];
    real<lower=0>  tau2[p];
  }
  model {
    sigma2 ~ inv_gamma(1,0.001);
    tau2int ~ uniform(0.01,1);
    alpha ~ normal(0, sqrt(100));
    for (k in 1:p) {
      gam2[k] ~ inv_gamma(0.5, 1);
      tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(sigma2*tau2[g[i]]));
    }
      for (i in 1:q) {
      betaint[i] ~ normal(0, sqrt(sigma2*(tau2[g1[i]]+tau2[g2[i]])/2*tau2int));
    }
    y ~ normal(alpha + X*beta + Xint*betaint,sqrt(sigma2));
  }
"

#Bayint* (with tau2int=1)
code_linint2_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> q;
    int<lower=0> n;
    matrix[n,p] X;
    matrix[n,q] Xint;
    real y[n];
    int g[p];
    int g1[q];
    int g2[q];
    real<lower=0.01, upper=1> tau2int;
  }
  parameters {
    real alpha;
    real<lower=0> sigma2;
    vector[p] beta;
    vector[q] betaint;
    real<lower=0>  gam2[p];
    real<lower=0>  tau2[p];
  }
  model {
    sigma2 ~ inv_gamma(1,0.001);
    alpha ~ normal(0, sqrt(100));
    for (k in 1:p) {
     gam2[k] ~ inv_gamma(0.5, 1);
      tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(sigma2*tau2[g[i]]));
    }
      for (i in 1:q) {
     betaint[i] ~ normal(0, sqrt(sigma2*sqrt(tau2[g1[i]]*tau2[g2[i]])*tau2int));
    }
    y ~ normal(alpha + X*beta + Xint*betaint,sqrt(sigma2));
  }
"

#Bayloc
code_lin_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> n;
    matrix[n,p] X;
    real y[n];
    int g[p];
  }
  parameters {
    real alpha;
    real<lower=0> sigma2;
    vector[p] beta;
    real<lower=0>  gam2[p];
    real<lower=0>  tau2[p];
  }
  model {
    sigma2 ~ inv_gamma(1,0.001);
    alpha ~ normal(0, sqrt(100));
    for (k in 1:p) {
      gam2[k] ~ inv_gamma(0.5, 1);
      tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(sigma2*tau2[g[i]]));
    }
    y ~ normal(alpha + X*beta,sqrt(sigma2));
  }
"

code_binint_brl_HalfCauchy <- "
  data { 
    int<lower=0> p;
    int<lower=0> q;
    int<lower=0> n;
    matrix[n,p] X;
    matrix[n,q] Xint;
    int<lower=0, upper=1> y[n];
    int g[p];
    int g1[q];
    int g2[q];
  }
  parameters {
    real alpha;
    real<lower=0.01, upper=1> tau2int;
    vector[p] beta;
    vector[q] betaint;
    real<lower=0>  gam2[p];
    real<lower=0>  tau2[p];
  }
  model {
    tau2int ~ uniform(0.01,1);
    alpha ~ normal(0, sqrt(100));
    for (k in 1:p) {
      gam2[k] ~ inv_gamma(0.5, 1);
      tau2[k] ~ inv_gamma(0.5, 1/gam2[k]);
    }
    for (i in 1:p) {
      beta[i] ~ normal(0, sqrt(tau2[g[i]]));
    }
      for (i in 1:q) {
      betaint[i] ~ normal(0, sqrt(sqrt(tau2[g1[i]]*tau2[g2[i]])*tau2int));
    }
    y ~ bernoulli_logit(alpha + X*beta + Xint*betaint);
  }
"


#make data chuncks
makeTraining <- function(data,chunksize,nfit=25,nrepeat=1){
  #chunksize <- 200;nfit=50;nrepeat<-1
  Y <- data[,1]
  X <- data[,-1] #remove response
  ssize <- length(Y)
  ssize
  ncov <- ncol(X)
  nset <- floor(ssize/chunksize)
  nset
  set.seed(34634)
  sets <- CVfolds(Y,model="linear", kfold=nset,nrepeat=nrepeat)
  if(nfit >= length(sets)) {print("increase nrepeat > nfit * chunksize/N");return(NULL)}
  
  #nfit <- 25
  set.seed(2152)
  nsetsfit <- sample(1:length(sets),nfit)
  
  
  trainind <- lapply(nsetsfit, function(ind){setind <- sets[[ind]]; return(setind)})
  return(trainind)
}

#map interactions to main effects
#cn <- colnames(mm)
#assumption: main effect names do not contain ':'
mapint <- function(cn){
  split <- sapply(cn,strsplit,split= ":")
  mains <- which(sapply(split,length)==1)
  splitint <- split[-mains]
  cnmain <- cn[mains]
  cnint <- cn[-mains]
  g1 <- unlist(lapply(splitint,function(spl){
    spl1 <- spl[1]
    return(match(spl1,cnmain))
  }))
  g2 <- unlist(lapply(splitint,function(spl){
    spl1 <- spl[2]
    return(match(spl1,cnmain))
  }))
  return(list(mains=mains,gs=rbind(g1,g2)))
}






