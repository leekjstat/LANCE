
rm(list=ls())

library(corrplot)
library(varband)
library(ggplot2)
library(reshape2)
library(gridExtra)

#######################################################################
# Main code
#######################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # setting the workspace to source file location
source("basic_fun.R")


#######################################################################
# 1-1. Simulation data
#######################################################################
p = 100; n = 100
A0min=0.4; A0max=0.6 # model parameters for generating the true Cholesky factor A0
D0min=2; D0max=5 # model parameters for generating the true diagonal matrix D0

# hyperparameters for LANCE prior
gamma = 0.1
alpha = 0.99
nu0 = 0
c1 = 1
c2 = 1.001


res.A0 = list()

res.LL = list()
res.YB = list()
res.LL.ROC = list()
res.YB.ROC = list()
res.LL.cv = list()
res.YB.cv = list()

n.sim = 10 # the number of repetition

for(sim in 1:n.sim){
  
  ########################################################################################################################################
  # 1. Data generation
  ########################################################################################################################################
  set.seed(sim)
  A0 = matrix(0, p,p)
  max.band = 0
  Rj = n/2 - 1
  for(i in 2:p){
    
    # # Model 1 : small band sizes
    # model.num = 1
    # max.band = 5
    # k0.i = min(sample(1:max.band, size=1), i-1)
    # vals = runif(n=k0.i, min=A0min, max=A0max) * sample(c(1,-1), size=k0.i, replace = TRUE)
    # A0[i, (i-k0.i):(i-1)] = vals
    
    # Model 2 in Yu and Bien : moderate band sizes
    model.num = 2
    max.band = 40
    if(runif(1) > 0.5 && (i-1)%%(p/5) > 0){
      k0.i = min(sample(1:((i-1)%%(p/5)), size=1), max.band)
      vals = runif(n=k0.i, min=A0min, max=A0max) * sample(c(1,-1), size=k0.i, replace = TRUE)
      A0[i, (i-k0.i):(i-1)] = vals
    }
    
    # # Model 3 (modified) : large band sizes
    # model.num = 3
    # max.band = 40
    # if(runif(1) > 0.5 && (i-1)%%(p/2) > 0){
    #   k0.i = min(sample(1:((i-1)%%(p/2)), size=1), max.band) # this part is modified! (to ensure k0.i < Rj = n/2)
    #   vals = runif(n=k0.i, min=A0min, max=A0max) * sample(c(1,-1), size=k0.i, replace = TRUE)
    #   A0[i, (i-k0.i):(i-1)] = vals
    # }
    
  }
  if(Rj < max.band) stop
  
  D0 = diag(runif(n = p, min = D0min, max = D0max))
  Omega0 = t(diag(p) - A0)%*%diag(1/diag(D0))%*%(diag(p) - A0)
  res.A0[[sim]] = A0
  
  # data generation (based on AR representation)
  X = matrix(0, nrow=n, ncol=p)
  X[,1] = rnorm(n, sd = sqrt(D0[1,1]))
  for(j in 2:p){
    mean.vec.j = X[, 1:(j-1)]%*%as.matrix(A0[j, 1:(j-1)])
    X[,j] = rnorm(n, mean = mean.vec.j, sd = sqrt(D0[j,j]))
  }
  
  
  ########################################################################################################################################
  # 2. Band test results
  ########################################################################################################################################
  
  #################################################################
  # CV-based inference 
  #################################################################
  # (1) Lee and Lin
  num = 5
  c2.vec = c(seq(from=-1.5, to=0.1, length.out=80), seq(from=0.13, to=1.5, length.out=10), seq(from=1.6, to=5, length.out=10))
  Rj = n/2 - 1
  
  res.LL.ROC[[sim]] = matrix(0, 2, length(c2.vec))
  dkj.hat.mat = matrix(0, p, min(Rj+1,j))
  for(jj in 1:length(c2.vec)){
    logpost.vec = list()
    for(j in 2:p){
      logpost.vec[[j]] = rep(0, min(Rj,j-1)+1)
      for(kj in 0:min(Rj,j-1)){
        if(jj == 1){
          res.jj = logpost.band(X, j, kj, gamma, nu0, c1, c2.vec[jj], alpha)
          logpost.vec[[j]][kj+1] = res.jj$val
          dkj.hat.mat[j, kj+1] = res.jj$dj
        }else{
          logpost.vec[[j]][kj+1] = - kj*log(1 + alpha/gamma)/2 - (alpha*n + nu0)*log(dkj.hat.mat[j, kj+1])/2 - kj*log(c1) - c2.vec[jj]*kj*log(p)
        }
      } # kj loop
    } # j loop
    MAP.bandwidths = rep(0,p)
    for(j in 2:p) MAP.bandwidths[j] = which.max(logpost.vec[[j]])-1
    
    MAP.A = matrix(0, p,p)
    for(j in 2:p){
      if(MAP.bandwidths[j] > 0) MAP.A[j, (j-MAP.bandwidths[j]):(j-1)] = 1
    }
    diag(MAP.A) = 0
    res.eval.LL = Evaluation((A0!=0)[lower.tri(A0)], MAP.A[lower.tri(A0)])
    res.LL.ROC[[sim]][1,jj] = res.eval.LL$Spe
    res.LL.ROC[[sim]][2,jj] = res.eval.LL$Sen
    
    # if(jj%%10 == 0)cat(jj%/%10 * 10, "th c2 is used. \n")
  } # jj loop
  res.LL[[sim]] = dkj.hat.mat
  
  res.LL.cv[[sim]] = logpost.band.cv(X, Rj, gamma, nu0, c1, c2.vec, alpha, num)
  
  
  #################################################################
  # (2) Yu & Bien (2017, JMLR)
  lambda.vec = seq(from=0.01, to=4, length.out = 100)
  res.YB.ROC[[sim]] = matrix(0, 2, length(lambda.vec))
  res.YB[[sim]] = list()
  for(jj in 1:length(lambda.vec)){
    res.YB.jj = varband(cov(X), lambda = lambda.vec[jj], diag(1/sqrt(diag(cov(X)))), w = FALSE)
    diag(res.YB.jj) = 0
    res.YB[[sim]][[jj]] =(abs(res.YB.jj) > 0.1^{10})
    
    res.eval = Evaluation((A0!=0)[lower.tri(A0)], res.YB[[sim]][[jj]][lower.tri(A0)])
    res.YB.ROC[[sim]][1,jj] = res.eval$Spe
    res.YB.ROC[[sim]][2,jj] = res.eval$Sen
    
    if(jj%%10 == 0)cat(jj%/%10 * 10, "th lambda is used. \n")
  }
  
  res.YB.cv[[sim]] = varband_cv(x = X, w = FALSE, lamlist = lambda.vec, nfolds = 5)
  
  
  #################################################################
  cat(sim, "th iteration has been completed. \n")
  #################################################################
}




main.ch = paste("Model",model.num," (n=",n,", p=",p,")", sep = "")
plot(res.YB.ROC[[1]][1,], res.YB.ROC[[1]][2,], xlim=c(0,1), ylim=c(0,1), xlab="specificity", ylab="sensitivity", 
     main=main.ch, lty=2, cex.main=2, cex.lab=1.5, type="l")
lines(res.LL.ROC[[1]][1,], res.LL.ROC[[1]][2,], col=2, type="l")

sim = 1
res = res.LL.cv[[sim]]$logpost.vec

MAP.bandwidths.cv = rep(0,p)
MAP.A.cv = matrix(0, p,p)
for(j in 2:p){
  MAP.bandwidths.cv[j] = which.max(res[[j]])-1
  if(MAP.bandwidths.cv[j] > 0) MAP.A.cv[j, (j-MAP.bandwidths.cv[j]):(j-1)] = 1
}
diag(MAP.A.cv) = 0
res1 = Evaluation((res.A0[[sim]]!=0)[lower.tri(res.A0[[sim]])], MAP.A.cv[lower.tri(res.A0[[sim]])])
LL.dot = c(res1$Spe, res1$Sen)

YB.cv = res.YB.cv[[sim]]$L_fit
diag(YB.cv) = 0
YB.cv =(abs(YB.cv) > 0.1^{10})
res2 = Evaluation((res.A0[[sim]]!=0)[lower.tri(res.A0[[sim]])], YB.cv[lower.tri(res.A0[[sim]])])
YB.dot = c(res2$Spe, res2$Sen)

points(LL.dot[1], LL.dot[2], col=2, cex=1.2, pch=16)
points(YB.dot[1], YB.dot[2], col=1, cex=1.2, pch=16)

for(sim in 2:n.sim){
  lines(res.YB.ROC[[sim]][1,], res.YB.ROC[[sim]][2,], lty=2, type="l")
  lines(res.LL.ROC[[sim]][1,], res.LL.ROC[[sim]][2,], col=2, type="l")
  
  res = res.LL.cv[[sim]]$logpost.vec
  MAP.bandwidths.cv = rep(0,p)
  MAP.A.cv = matrix(0, p,p)
  for(j in 2:p){
    MAP.bandwidths.cv[j] = which.max(res[[j]])-1
    if(MAP.bandwidths.cv[j] > 0) MAP.A.cv[j, (j-MAP.bandwidths.cv[j]):(j-1)] = 1
  }
  diag(MAP.A.cv) = 0
  res1 = Evaluation((res.A0[[sim]]!=0)[lower.tri(res.A0[[sim]])], MAP.A.cv[lower.tri(res.A0[[sim]])])
  LL.dot = c(res1$Spe, res1$Sen)
  
  YB.cv = res.YB.cv[[sim]]$L_fit
  diag(YB.cv) = 0
  YB.cv =(abs(YB.cv) > 0.1^{10})
  res2 = Evaluation((res.A0[[sim]]!=0)[lower.tri(res.A0[[sim]])], YB.cv[lower.tri(res.A0[[sim]])])
  YB.dot = c(res2$Spe, res2$Sen)
  
  points(LL.dot[1], LL.dot[2], col=2, cex=1.2, pch=16)
  points(YB.dot[1], YB.dot[2], col=1, cex=1.2, pch=16)
}
legend("bottomleft", legend=c("LL", "LL (cv)", "YB", "YB (cv)"),
       bty="n", # remove border
       col=c(2,2,1,1),
       lty=c(1,1,2,2), # line type
       pch=c(NA,16,NA,16),# marker shape
       cex=1.6, # shape size
       lwd=c(2,0,2,0), # line width
       text.font = 2) # style of the font




