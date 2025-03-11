
library(readr)
library(lavaan)
library(lessSEM)
library(glmnet)

EZest <- function(A, A0, B, B0, Gamma, Phi1, Phi2, Phi3, X, Y){
  q <- ncol(X)
  p <- ncol(Y)
  s <- ncol(B)
  m <- ncol(A)
  n <- nrow(X)
  
  sigma21 <-  rbind( matrix( cbind( t(B), t(Gamma)%*% t(A)), nrow = s, ncol = q+p)    , matrix(cbind(Gamma%*%t(B),(Phi3 + Gamma%*%t(Gamma))%*%t(A)), nrow = m, ncol = q+p) )  ## m (m for Z)  s (s for W): Z = Gamma W + E
  
  sigma11 <- matrix(cbind(rbind((Phi1 + B%*%t(B)), A%*%Gamma%*%t(B)), rbind(B%*%t(A%*%Gamma),
                                                                                        (Phi2 + A%*%(Phi3 + Gamma%*%t(Gamma))%*%t(A)))), nrow = q+p, ncol = q+p)
  
  sigma22 <-   rbind( matrix( cbind( diag(s), t(Gamma) ), nrow = s, ncol = (m + s)   )  ,matrix( cbind( Gamma  , Phi3 + Gamma%*%t(Gamma) ), nrow = m, ncol = (m + s)   ) ) ## (m +s) x (m +s )
  
  condvarWZ <- sigma22 - sigma21%*%solve(sigma11)%*%t(sigma21)  ##  (m +s) x (m +s)
  condvar  <- condvarWZ[s+(1:m), s+(1:m) ]  ## m x m    
  EWZ <- sigma21%*%solve(sigma11)%*%(rbind(t(X),t(Y))-matrix(rbind(B0,A0), nrow = q+p, ncol = n, byrow = F)) ## (m + s) x n
  EZ <-  matrix(EWZ[s+(1:m), ], nrow = m)  ## p x n
  EZ 
}

COVIDiSTRESS_global_survey_May_30_2020_final_cleaned_file_ <- read_csv("COVIDiSTRESS global survey May 30 2020 (___final cleaned file___).csv")

BICs <- c()

Country.PSEs <- data.frame("Country" = as.character(), "PSE" = as.numeric(), "Opt_s" = as.numeric(), "Opt_m" = as.numeric())
Country.mlm.PSEs <- data.frame("Country" = as.character(), "PSE" = as.numeric())
Country.lasso.PSEs <- data.frame("Country" = as.character(), "PSE" = as.numeric())
Country.lav.PSEs <- data.frame("Country" = as.character(), "PSE" = as.numeric())
Country.less.PSEs <- data.frame("Country" = as.character(), "PSE" = as.numeric())
Country.sparsity <- data.frame("Country" = as.character(), "prop. B" = as.numeric(), "prop.A" = as.numeric())

Countries <- c("Australia", "Austria", "Brazil", "Canada", "Croatia", "Germany", "Greece",
               "Ireland", "Indonesia", "Malaysia", "Netherlands", "Pakistan", "Panama", "Philippines",
               "Poland", "Portugal", "Serbia", "Slovakia", "South Korea", "Spain",
               "Switzerland", "Turkey", "United Kingdom", "United States")

for(Country in Countries){
  
  Output <- readRDS(paste0(Country, ".Results.rds"))
  
  set.seed(2435)
  
  if(Country == "South Korea"){
    Country.dat <- COVIDiSTRESS_global_survey_May_30_2020_final_cleaned_file_[COVIDiSTRESS_global_survey_May_30_2020_final_cleaned_file_$Country == "Korea, South",]
  }else{
    Country.dat <- COVIDiSTRESS_global_survey_May_30_2020_final_cleaned_file_[COVIDiSTRESS_global_survey_May_30_2020_final_cleaned_file_$Country == Country,]
  }
  Country.dat <- Country.dat[!is.na(Country.dat$Country),]
  Country.dat <- Country.dat[,c(22:48,49:54,121:136)]
  Country.dat <- Country.dat[complete.cases(Country.dat),]
  X.Country <- Country.dat[,c(1:27)]
  Y.Country <- Country.dat[,c(28:49)]
  Country.ind <- sample(seq_len(nrow(X.Country)), size = floor(0.60*nrow(X.Country)))
  X.Country.train <- X.Country[Country.ind,]
  X.Country.test <- X.Country[-Country.ind,]
  Y.Country.train <- Y.Country[Country.ind,]
  Y.Country.test <- Y.Country[-Country.ind,]
  EZ.test <- EZest(Output$A_estimates$A_final_pars$A, Output$A_estimates$A_final_pars$A0, Output$B_estimates$B_final_pars$B, Output$B_estimates$B_final_pars$B0,
                   Output$A_estimates$A_final_pars$Gamma, Output$B_estimates$B_final_pars$Phi1, Output$A_estimates$A_final_pars$Phi2, Output$A_estimates$A_final_pars$Phi3, X.Country.test, Y.Country.test)

  Country.PSEs[which(Countries == Country),1] <- Country
  Country.PSEs[which(Countries == Country),2] <- 1/(ncol(Y.Country.test)*nrow(Y.Country.test))*sum((Y.Country.test - t(matrix(Output$A_estimates$A_final_pars$A0,nrow = ncol(Y.Country.test), ncol = nrow(Y.Country.test), byrow = F) + Output$A_estimates$A_final_pars$A%*%EZ.test))^2)
  Country.PSEs[which(Countries == Country),3] <- Output$B_estimates$s_opt
  Country.PSEs[which(Countries == Country),4] <- Output$A_estimates$m_opt

Country.lm <- lm(as.matrix(Y.Country.train) ~ as.matrix(X.Country.train))
MLM.pred <- cbind(rep(1, nrow(X.Country.test)),as.matrix(X.Country.test))%*%Country.lm$coefficients
Country.mlm.PSEs[which(Countries == Country),1] <- Country
Country.mlm.PSEs[which(Countries == Country),2] <- 1/(ncol(Y.Country.test)*nrow(Y.Country.test))*sum((MLM.pred - Y.Country.test)^2)

Country.lasso.cv <- cv.glmnet(as.matrix(X.Country.train), as.matrix(Y.Country.train), family = "mgaussian")
Lasso.pred <- predict(Country.lasso.cv, newx = as.matrix(X.Country.test), s = "lambda.min")
Country.lasso.PSEs[which(Countries == Country),1] <- Country
Country.lasso.PSEs[which(Countries == Country),2] <- 1/(ncol(Y.Country.test)*nrow(Y.Country.test))*sum((Lasso.pred - Y.Country.test)^2)

Lav.output <- readRDS(paste0(Country, ".lavaan.Results.rds"))
Lav.pred <- lavPredictY(Lav.output, cbind(X.Country.test, Y.Country.test), ynames = colnames(Y.Country.test), xnames = colnames(X.Country.test))
Country.lav.PSEs[which(Countries == Country),1] <- Country
Country.lav.PSEs[which(Countries == Country),2] <- 1/(ncol(Y.Country.test)*nrow(Y.Country.test))*sum((Lav.pred - Y.Country.test)^2)

if(file.exists(paste0(Country, ".lessSEM.Results.rds"))){
Less.output <- readRDS(paste0(Country, ".lessSEM.Results.rds"))
Less.pred <- lavPredictY(lessSEM2Lavaan(Less.output, criterion = "BIC"), cbind(X.Country.test, Y.Country.test), ynames = colnames(Y.Country.test), xnames = colnames(X.Country.test))
Country.less.PSEs[which(Countries == Country),1] <- Country
Country.less.PSEs[which(Countries == Country),2] <- 1/(ncol(Y.Country.test)*nrow(Y.Country.test))*sum((Less.pred - Y.Country.test)^2)
}

nzeroB <- sum(Output$B_estimates$B_final_pars$B == 0) - (ncol(Output$B_estimates$B_final_pars$B)*(ncol(Output$B_estimates$B_final_pars$B)-1))/2
nfreeB <- ncol(Output$B_estimates$B_final_pars$B)*(2*nrow(Output$B_estimates$B_final_pars$B) - ncol(Output$B_estimates$B_final_pars$B) - 1)/2
nzeroA <- sum(Output$A_estimates$A_final_pars$A == 0) - (ncol(Output$A_estimates$A_final_pars$A)*(ncol(Output$A_estimates$A_final_pars$A)-1))/2
nfreeA <- ncol(Output$A_estimates$A_final_pars$A)*(2*nrow(Output$A_estimates$A_final_pars$A) - ncol(Output$A_estimates$A_final_pars$A) - 1)/2
Country.sparsity[which(Countries == Country),1] <- Country
Country.sparsity[which(Countries == Country),2] <- nzeroB/(nfreeB)
Country.sparsity[which(Countries == Country),3] <- nzeroA/(nfreeA)

}
