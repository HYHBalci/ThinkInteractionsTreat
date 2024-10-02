

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

#auxiliary function to retrieve cred int from posterior
postsum <- function(posters) {
  qs <- quantile(posters,p=c(0.025,0.975))
  mn <- mean(posters)
  return(c(qs[1],mn,qs[2]))
}

#posterior of age + age*gender
post2f <- function(j,k,data) {
  #k<-"bmi";j<-"age"; data <- testdatx
  cn <- colnames(data)
  int2main <- mapint(cn)
  mains <- int2main$main
  cnmains <- cn[mains]
  ints <- int2main$gs
  if(is.character(j)) {
    if(is.na(match(j,cnmains)) | is.na(match(k,cnmains))) return("Error: use correct main effect name")
    else {
      j <- which(cn == j); 
      if(is.character(k)) {
        whk <- which(cn == k); 
        k <- which(apply(ints,2,function(vec) max(sum(vec==c(j,whk)),sum(vec==c(whk,j))))==2)
      }
    }
  }
  ints <- int2main$gs
  whints <- ints[,k]
  whfocus <- whints[whints==j]
  if(length(whfocus)==0) return("Error: interaction does not match with main effect")
  else {
    whother <- setdiff(whints,whfocus)
    lev <- unique(data[,whother])
    postjk <- sapply(lev,function(xk) return(betapostb[,j] + betaposti[,k]*xk))
    posterint <- rbind(lev,apply(postjk,2,postsum))
    return(posterint)
  }
}

#personalized unit change posterior
postallf <- function(j,data,betapostb,betaposti) {
  #j<-"age"; data <- Xtest
  cn <- colnames(data)
  int2main <- mapint(cn)
  mains <- int2main$main
  cnmains <- cn[mains]
  ints <- int2main$gs
  if(is.character(j)) {
    if(is.na(match(j,cnmains))) return("Error: use correct main effect name")
    else {
      j <- which(cn == j)
    }
  }
  ints <- int2main$gs
  whints <- which(apply(ints,2,function(ab) is.element(j,ab)))
  whothers <- apply(ints[,whints],2,function(ab) setdiff(ab,j))
  postjk <- betapostb[,j] + t(as.matrix(data[,whothers]) %*% t(betaposti[,whints]))
  posterint <- apply(postjk,2,postsum)
  return(posterint)
}

#shapley values and their uncertainty
shapley <- function(js,data,exx,ex,betapostb,betaposti) {
data <- as.matrix(data)
cn <- colnames(data)
int2main <- mapint(cn)
mains <- int2main$main
cnmains <- cn[mains]
ints <- int2main$gs
if(is.character(js)) {
  if(is.na(match(js,cnmains))) return("Error: use correct main effect name")
  else {
    js <- which(cn == js)
  }
}
ints <- int2main$gs
shapleylinj <- 0;  shapleyintj <- 0;
for(j in js){
  whints <- which(apply(ints,2,function(ab) is.element(j,ab)))
  whothers <- apply(ints[,whints],2,function(ab) setdiff(ab,j))
  shapleylinj <- shapleylinj + (data[,j,drop=FALSE]-ex[j]) %*% t(betapostb[,j,drop=FALSE])
  shapleyintj <- shapleyintj +1/2*( (data[,j] - ex[j])*data[,whothers,drop=FALSE] + 
                                      t(t(data[,j,drop=FALSE] %*% matrix(ex[whothers],nrow=1)) -exx[whothers,j])) %*% t(betaposti[,whints])
}
shlinpost  <- t(apply(shapleylinj,1,postsum))
shintpost <- t(apply(shapleyintj,1,postsum))
shtotpost <- t(apply(shapleylinj + shapleyintj,1,postsum))
ress <- list(tot=shtotpost,lin=shlinpost,int=shintpost)
return(ress)
}

shapleys <- function(trainmain,testdat, betapostb,betaposti, grouped = list()){
  #grouped <- list(etn=3:6);trainmain <- traindatxmain;testdat <- testdatx
  EXX <- cov(trainmain)
  EX <- apply(trainmain,2,mean)
  nmain <- ncol(betapostb)
  cn <- colnames(trainmain)
  if(length(grouped)==0) {groups = as.list(1:nmain); names(groups) <- cn}  else {
    gr <- unlist(grouped)
    notin <- setdiff(1:nmain,gr)
    nin <- cn[notin]
    group1 <- as.list(notin)
    names(group1) <- nin
    groups <- c(group1, grouped)
    od <- order(unlist(lapply(groups,min)))
    groups <- groups[od]
  }
  res <- lapply(groups,shapley,data=testdat,exx=EXX,ex=EX,betapostb=betapostb,betaposti=betaposti)
  return(res)
}







