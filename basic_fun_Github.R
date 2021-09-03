

#######################################################################
# log-posterior probability for varying bandwidths Lee and Lin
#######################################################################
logpost.band <- function(X, j, kj, gamma, nu0, c1, c2, alpha){
   # X: data matrix
   # j: jth colunm
   # kj: bandwidth 
   # gamma, nu0, c1, c2: hyperparameters
   # alpha: alpha-fractional parameter
   
   n = nrow(X)
   p = ncol(X)
   
   dkj.hat = dhat(j, kj, X)
   logpost.val = - kj*log(1 + alpha/gamma)/2 - (alpha*n + nu0)*log(dkj.hat)/2 - kj*log(c1) - c2*kj*log(p)
   return( list(val = logpost.val, dj = dkj.hat) )
}

logpost.band.cv <- function(X, Rj, gamma, nu0, c1, c2.vec, alpha, num){
   # Bayesian CV
   # X: data matrix
   # Rj: maximum bandwidth 
   # gamma, nu0, c1, c2: hyperparameters
   # alpha: alpha-fractional parameter
   # num: number of random partitions used for CV
   
   n = nrow(X)
   p = ncol(X)
   c2.len = length(c2.vec)
   R.vec = rep(0, c2.len)
   
   
   for(ind in 1:num){    
      i.ind = sample(1:n, size=n/2)
      X.train = X[i.ind,]
      X.test = X[-i.ind,]
      
      post.vec = list()
      logpost.vec = list()
      dkj.i.vec = list()
      dkj.hat.mat = list()
      
      for(jj in 1:length(c2.vec)){
         
         log.pred = 0
         for(j in 2:p){
            post.vec[[j]] = rep(0, min(Rj,j-1)+1)
            logpost.vec[[j]] = rep(0, min(Rj,j-1)+1)
            
            if(jj == 1){
               dkj.i.vec[[j]] = rep(0, min(Rj,j-1)+1)
               dkj.hat.mat[[j]] = rep(0, min(Rj,j-1)+1)
            }
            
            # calculating (1) \pi_\alpha(kj | X.train) = post.vec[[j]][kj+1] and (2) dhat(j, kj, X.test) = dkj.i.vec[[j]][kj+1] for given c2.vec[[jj]]
            for(kj in 0:min(Rj,j-1)){
               if(jj == 1){
                  res.jj = logpost.band(X.train, j, kj, gamma, nu0, c1, c2.vec[jj], alpha)
                  logpost.vec[[j]][kj+1] = res.jj$val
                  dkj.hat.mat[[j]][kj+1] = res.jj$dj
                  dkj.i.vec[[j]][kj+1] = dhat(j, kj, X.test)
                  # cat(kj, res.jj$dj,"\n")
               }else{
                  logpost.vec[[j]][kj+1] = - kj*log(1 + alpha/gamma)/2 - (alpha*nrow(X.train) + nu0)*log(dkj.hat.mat[[j]][kj+1])/2 - kj*log(c1) - c2.vec[jj]*kj*log(p)
               }
               
            } # kj loop
            
            post.vec[[j]] = exp(logpost.vec[[j]] - max(logpost.vec[[j]]))
            post.vec[[j]] = post.vec[[j]]/sum(post.vec[[j]]) # normalizing the posterior prob.
            
            # log of predictive density of Xj in X.test
            log.pred = log.pred + log(sum( (1+1/gamma)^{-(0:min(Rj,j-1))/2} * dkj.i.vec[[j]]^{-(1*nrow(X.test) + nu0)/2} * post.vec[[j]], na.rm = TRUE ))
         } # j loop
         
         R.vec[jj] = R.vec[jj] + log.pred
         # cat(jj,"th c2, and", ind,"th ind. \n")
      } # jj loop
      
      # if(ind%%10 == 0)cat(ind%/%10 * 10, "th iteration has completed. \n")
      cat(ind, "th random split .. \n")
   } # ind loop
   
   
   
   max.ind = which.max(R.vec)
   c2.sel = c2.vec[max.ind]
   
   logpost.vec = list()
   dkj.hat.mat = list()
   for(j in 2:p){
    logpost.vec[[j]] = rep(0, min(Rj,j-1)+1)
    dkj.hat.mat[[j]] = rep(0, min(Rj,j-1)+1)
      for(kj in 0:min(Rj,j-1)){
         res.jj = logpost.band(X, j, kj, gamma, nu0, c1, c2.sel[1], alpha)
         logpost.vec[[j]][kj+1] = res.jj$val
         dkj.hat.mat[[j]][kj+1] = res.jj$dj
      }
   }
   
   return( list(logpost.vec = logpost.vec, dkj.hat.mat = dkj.hat.mat, c2.sel = c2.sel, R.vec = R.vec) )
}

#######################################################################
# Auxiliary functions 
#######################################################################

# Z_j(k) : n X k matrix
Z <- function(j, k, X){
   ind = c(max(1, j-k):(j-1))
   res = X[, ind]
   if(j == 2 || k == 1) res = matrix(res)
   return(res)
}

# \hat{Var}(Z_j(k)) : k X k matrix
VhatZ <- function(j, k, X){
   Zj = Z(j, k, X)  
   res = t(Zj)%*%Zj/nrow(X)
   return(res)
}

# \hat{d}_{jk} 
dhat <- function(j, k, X){
   n = nrow(X)
   Xj = X[, j]
   if(j == 1 || k == 0){
      res = sum(Xj^2)/n   # scalar
   }else{
      Zj = Z(j, k, X)
      VhatX = sum(Xj^2)/n   # scalar
      Covhat = t(Zj)%*%matrix(Xj)/n   # k X 1 matrix
      
      Vhat = VhatZ(j, k, X)
      L = eigen(Vhat)$vec
      D = eigen(Vhat)$val; nd = length(D)
      D = diag(x=D, nrow=nd)
      Vhat.inv = L%*%diag(1/diag(D), nrow=nd)%*%t(L)
      
      res = VhatX - t(Covhat)%*% Vhat.inv %*% Covhat
      if(res <= 0){
         ahat = Vhat.inv %*% Covhat
         
         # res = t(as.matrix(c(-ahat, 1))) %*% VhatZ(j+1,k+1,X) %*% as.matrix(c(-ahat, 1))
         ahat1 = as.matrix(c(-ahat, 1))
         Zj1 = Z(j+1, k+1, X)/nrow(X)
         res = sum((Zj1 %*% ahat1)^2)
         
      }
   }
   return(res)
}


Evaluation <- function(beta1, beta2){
   true.index <- which(beta1==1)
   false.index <- which(beta1==0)
   positive.index <- which(beta2==1)
   negative.index <- which(beta2==0)
   
   TP <- length(intersect(true.index,positive.index))
   FP <- length(intersect(false.index,positive.index))
   FN <- length(intersect(true.index,negative.index))
   TN <- length(intersect(false.index,negative.index))
   
   
   Precision <- TP/(TP+FP)
   if((TP+FP)==0) Precision <- 1
   Recall <- TP/(TP+FN)
   if((TP+FN)==0) Recall <- 1
   Sensitivity <- Recall
   Specific <- TN/(TN+FP)
   if((TN+FP)==0) Specific <- 1
   
   if(TP*TN - FP*FN == 0){
   	MCC = NA
   }else if(TP*TN - FP*FN > 0){
   	logMCC = log(TP*TN - FP*FN) - 0.5 * (log(TP+FP) + log(TP+FN) + log(TN+FP) + log(TN+FN))
   	MCC = exp(logMCC)
   }else{
   	logMCC = log(-TP*TN + FP*FN) - 0.5 * (log(TP+FP) + log(TP+FN) + log(TN+FP) + log(TN+FN))
   	MCC = -exp(logMCC)
   }
   
   
   return(list(Precision=Precision,Sensitivity=Sensitivity,Specific=Specific, MCC=MCC, TP=TP,FP=FP,TN=TN,FN=FN))
}


predict.err <- function(est.L, X.test){
   p = ncol(X.test)
   n = nrow(X.test)
   
   res = rep(0, n)
   for(i in 1:n){
      x.i = matrix(X.test[i,])
      res[i] = sum( (est.L %*% x.i)^2 )/(p-1)
   }
   
   return(res)
}





