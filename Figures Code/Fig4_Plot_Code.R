library(readr)
library(gdata)
library(glmnet)
library(lavaan)
library(lessSEM)
library(ggplot2)

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

ResTable <- replicate(6,data.frame("n" = double(), "PSE" = as.double(), "Model" = character()),simplify=F)

set.seed(2435)

source(file = "Model1_Meta_Param.R")
source(file = "Model1_Matrix_Param.R")
source(file = "Data_Generation.R")

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("Lavrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])
    
    datset <- as.numeric(gsub(".*rep(.+)n.*", "\\1", datfiles[j]))
    set.seed(datset)
    
    dat <- Data_Generation(n = i, Model_Params=Model_Params)
    
    Split.ind <- sample(seq_len(nrow(dat$X)), size = floor(0.60*nrow(dat$X)))
    X.train <- dat$X[Split.ind,]
    X.test <- dat$X[-Split.ind,]
    Y.train <- dat$Y[Split.ind,]
    Y.test <- dat$Y[-Split.ind,]
    
    colnames(X.test) <- sprintf("X%d", 1:q)
    colnames(Y.test) <- sprintf("Y%d", 1:p)
    
    ResTable[[1]][counter, 1] <- i
    ResTable[[1]][counter, 2] <- 1/(p*nrow(Y.test))*sum((lavPredictY(Output, cbind(X.test, Y.test), ynames = colnames(Y.test), xnames = colnames(X.test)) - Y.test)^2)
    ResTable[[1]][counter, 3] <- "Lavaan"
    counter <- counter + 1
  }
}

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("Overallrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])
    
    datset <- as.numeric(gsub(".*rep(.+)n.*", "\\1", datfiles[j]))
    set.seed(datset)
    
    dat <- Data_Generation(n = i, Model_Params=Model_Params)
    
    Split.ind <- sample(seq_len(nrow(dat$X)), size = floor(0.60*nrow(dat$X)))
    X.train <- dat$X[Split.ind,]
    X.test <- dat$X[-Split.ind,]
    Y.train <- dat$Y[Split.ind,]
    Y.test <- dat$Y[-Split.ind,]
    
    EZ.test <- EZest(Output$A_estimates$A_final_pars$A, Output$A_estimates$A_final_pars$A0, Output$B_estimates$B_final_pars$B, Output$B_estimates$B_final_pars$B0,
                     Output$A_estimates$A_final_pars$Gamma, Output$B_estimates$B_final_pars$Phi1, Output$A_estimates$A_final_pars$Phi2,Output$A_estimates$A_final_pars$Phi3, X.test, Y.test)
    
    ResTable[[2]][counter, 1] <- i
    ResTable[[2]][counter, 2] <- 1/(p*nrow(Y.test))*sum((Y.test - t(matrix(Output$A_estimates$A_final_pars$A0,nrow = ncol(Y.test), ncol = nrow(Y.test), byrow = F) + Output$A_estimates$A_final_pars$A%*%EZ.test))^2)
    ResTable[[2]][counter, 3] <- "JSFRM"
    
    MLM.lm.fit <- cbind(rep(1, nrow(X.test)), X.test)%*%Output[[2]]$coefficients
    ResTable[[3]][counter, 1] <- i
    ResTable[[3]][counter, 2] <- 1/(p*nrow(Y.test))*sum((MLM.lm.fit - Y.test)^2)
    ResTable[[3]][counter, 3] <- "MLM"
    
    Lasso.fit <- predict(Output[[3]], newx = X.test, s = "lambda.min")[,,1]
    ResTable[[4]][counter, 1] <- i
    ResTable[[4]][counter, 2] <- 1/(p*nrow(Y.test))*sum((Lasso.fit - Y.test)^2)
    ResTable[[4]][counter, 3] <- "Lasso"
    
    B.comp <- Output[[4]][[1]]
    A.comp <- Output[[4]][[2]]
    G.comp <- Output[[4]][[3]]
    EW.test <- EWest(B.comp$Bopt$B, B.comp$Bopt$B0, B.comp$Bopt$Phi1, X.test)
    
    ResTable[[5]][counter, 1] <- i
    ResTable[[5]][counter, 2] <- 1/(p*nrow(Y.test))*sum((t(matrix(A.comp$Bopt$B0, nrow = ncol(Y.test), ncol = nrow(Y.test), byrow = F) + A.comp$Bopt$B%*%(t(G.comp[-1,])%*%EW.test)) - Y.test)^2)
    ResTable[[5]][counter, 3] <- "IFR"
    
    counter <- counter + 1
  }
}

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("Lessrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])
    Output <- lessSEM2Lavaan(Output, criterion = "BIC")
    
    datset <- as.numeric(gsub(".*rep(.+)n.*", "\\1", datfiles[j]))
    set.seed(datset)
    
    dat <- datgen(params1Astotal, params1Bstotal, params1G, i, varvecB, varvecA, varvecG)
    
    Split.ind <- sample(seq_len(nrow(dat$X)), size = floor(0.60*nrow(dat$X)))
    X.train <- dat$X[Split.ind,]
    X.test <- dat$X[-Split.ind,]
    Y.train <- dat$Y[Split.ind,]
    Y.test <- dat$Y[-Split.ind,]
    
    colnames(X.test) <- sprintf("X%d", 1:q)
    colnames(Y.test) <- sprintf("Y%d", 1:p)
    
    ResTable[[6]][counter, 1] <- i
    ResTable[[6]][counter, 2] <- 1/(p*nrow(Y.test))*sum((lavPredictY(Output, cbind(X.test, Y.test), ynames = colnames(Y.test), xnames = colnames(X.test)) - Y.test)^2)
    ResTable[[6]][counter, 3] <- "LessSEM"
    counter <- counter + 1
  }
}

ResTable <- do.call(rbind, ResTable)

ggplot(data = ResTable, aes(x = as.factor(n), y = PSE, fill = Model)) + geom_boxplot(outlier.shape = NA) + xlab("Sample Size") + ylab("PSE")  + theme_bw() + scale_fill_grey() + scale_fill_discrete(labels=c("IFR", "JSFRM", "Lasso", "Lavaan", "LessSEM", "MLM"))
