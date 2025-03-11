args <- commandArgs(TRUE)
datset <- as.numeric(args[1]) 
sampsize <- as.numeric(args[2])

library(glmnet)
library(parallel)

EWest <- function(B, B0, Phi1, X){
  q <- ncol(X)
  xk_sig <- ncol(B)
  n <- nrow(X)
  
  invmodv <-Matrix::solve(Phi1 + tcrossprod(B))
  condvar <- diag(xk_sig) - t(B)%*%invmodv%*%B
  
  EW <- t(B)%*%invmodv%*%(t(X)-matrix(B0, nrow = q, ncol = n, byrow = F))
  EW
}

EZest <- function(A, A0, B, B0, Gamma, Phi1, Phi2, Phi3, X, Y){
  q <- ncol(X)
  p <- ncol(Y)
  xk_sig <- ncol(B)
  k_sig <- ncol(A)
  n <- nrow(X)
  
  sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = xk_sig, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = k_sig, ncol = q+p) )  ## k_sig (m for Z)  xk_sig (s for W): Z = Gamma W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + tcrossprod(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                 (Phi2 + A%*%(Phi3 + tcrossprod(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( diag(xk_sig), t(Gamma) ), nrow = xk_sig, ncol = (k_sig + xk_sig)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = k_sig, ncol = (k_sig + xk_sig)   ) ) ## (m +s) x (m +s )
  
  condvarWZ <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)  ##  (m +s) x (m +s)
  condvar  <- condvarWZ[xk_sig+(1:k_sig), xk_sig+(1:k_sig) ]  ## m x m    
  EWZ <- sigma21%*%solve(sigma11)%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
  EZ <-  matrix(EWZ[xk_sig+(1:k_sig), ], nrow = k_sig)  ## p x n
  EZ 
}

set.seed(2435)

source(file = "Model1_Meta_Param.R")
source(file = "Model1_Matrix_Param.R")

source(file = "Data_Generation.R")

set.seed(datset)

dat <- Data_Generation(n = sampsize, Model_Params=Model_Params)

source(file = "B_inits.R")
source(file = "Estep_X.R")
source(file = "Mstep_X.R")
source(file = "logLik_X.R")
source(file = "Convergence_check.R")
source(file = "BIC_X.R")
source(file = "EMAlgBAdLassoCV.R")
source(file = "OverallBAlg.R")

Singular_ErrorX <- function(tuningpB, Data, initial_Model_Params, weights){
  return(tryCatch(EMAlgBAdLassoCV(tuningpB, Data, initial_Model_Params, weights), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "EMAlgBAdLassoCV"))

Best <- OverallBAlg(MyData$X, s_seq = 3)

stopCluster(cl)

source(file = "A_inits.R")
source(file = "Estep_Y.R")
source(file = "Mstep_Y.R")
source(file = "logLik_Y.R")
source(file = "BIC_Y.R")
source(file = "EMAlgAGammaAdLassoCV.R")
source(file = "OverallAGAlg.R")

Singular_ErrorY <- function(tuningpA, Data, Best, initial_Model_Params, weights){
  return(tryCatch(EMAlgAGammaAdLassoCV(tuningpA, Data, Best, initial_Model_Params, weights), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_Y", "Mstep_Y", "logLik_Y", "Convergence_check", "BIC_Y", "EMAlgAGammaAdLassoCV"))

Aest <- OverallAGAlg(MyData, Best, m_seq = 2)

stopCluster(cl)

Results.EZ <- EZest(Output$A_estimates$A_final_pars$A, Output$A_estimates$A_final_pars$A0, Output$B_estimates$B_final_pars$B, Output$B_estimates$B_final_pars$B0,
                    Output$A_estimates$A_final_pars$Gamma, Output$B_estimates$B_final_pars$Phi1, Output$A_estimates$A_final_pars$Phi2,Output$A_estimates$A_final_pars$Phi3, dat$X, dat$Y)

Results.PSE <- sum((1/(nrow(dat$Y)*ncol(dat$Y)))*(dat$Y - t(matrix(Results$A_estimates$A_final_pars$A0,nrow = ncol(dat$Y), ncol = nrow(dat$Y), byrow = F) + Results$A_estimates$A_final_pars$A%*%Results.EZ))^2)

lin.mod <- lm(dat$Y ~ dat$X)
lm.fit <- cbind(rep(1, nrow(dat$X)),dat$X)%*%lin.mod$coefficients
lm.PSE <- sum((1/(nrow(dat$Y)*ncol(dat$Y)))*(dat$Y - lm.fit)^2)

lasso.mod <- cv.glmnet(dat$X, dat$Y, family = "mgaussian")
lasso.fit <- predict(lasso.mod, newx = dat$X, s = "lambda.min")
lasso.PSE <- sum((1/(nrow(dat$Y)*ncol(dat$Y)))*(dat$Y - lasso.fit[,,1])^2)

B.comp <- OverallBAlg(dat$X, s_seq = 3)
A.comp <- OverallBAlg(dat$Y, s_seq = 2)

B.comp.EW <- EWest(B.comp$B_final_pars$B, B.comp$B_final_pars$B0, B.comp$B_final_pars$Phi1, dat$X)
A.comp.EW <- EWest(A.comp$B_final_pars$B, A.comp$B_final_pars$B0, A.comp$B_final_pars$Phi1, dat$Y)

G.comp <- lm(t(A.comp.EW) ~ t(B.comp.EW))$coefficients

XY.factor.PSE <- sum((1/(nrow(dat$Y)*ncol(dat$Y)))*(dat$Y - t(matrix(A.comp$B_final_pars$B0,nrow = ncol(dat$Y), ncol = nrow(dat$Y), byrow = F) + A.comp$B_final_pars$B%*%(t(G.comp[-1,])%*%B.comp.EW)))^2)

Output <- list(Results.PSE, lm.PSE, lasso.PSE, XY.factor.PSE)

save(Output, file=paste("Overallrep", datset, "n", sampsize, ".Rdata", sep=""))
