#Model4_Table_Code.R

ResTable <- replicate(11,data.frame(matrix(0,2000,2), repnames = character(2000)),simplify=F)
count <- 0
iter <- 0
optlatB <- NULL
optlatA <- NULL
opttunB <- NULL
opttunA <- NULL
AMSElist <- NULL
BMSElist <- NULL
GMSElist <- NULL
TPRA <- NULL
TNRA <- NULL
TPRB <- NULL
TNRB <- NULL
itb1 <- rep(0, 4)
ita1 <- rep(0, 4)
itg1 <- rep(0, 4)
itb2 <- rep(0, 2)
ita2 <- rep(0, 2)
itg2 <- rep(0, 2)

set.seed(2435)

source(file = "Model4_Meta_Param.R")
source(file = "Model45_Matrix_Param.R")

totreps <- rep(0, 4)

for(i in c(50, 60, 70, 80)){
  datfiles <- list.files(pattern=paste0("Model4rep.*n",i,".rds"))
  totreps[which(i == c(50, 60, 70, 80))] <- length(datfiles)
  if(length(datfiles) > 0){
    for(k in 1:length(datfiles)){
      opttunB <- NULL
      opttunA <- NULL
      BMSElist <- NA
      AMSElist <- NA
      GMSElist <- NA
      TNRA <- NA
      TPRA <- NA
      TNRB <- NA
      TPRB <- NA
      Results <- readRDS(datfiles[k])
      optlatB <- Results$B_estimates$s_opt
      optlatA <- Results$A_estimates$m_opt
      opttunB <- Results$B_estimates$tuningpB
      opttunA <- Results$A_estimates$tuningpA
      if(optlatB == meta_param$s){
        TNRB <- 0
        TPRB <- 0
        for(u in 1:meta_param$s){
          for(v in 1:meta_param$s){
            if(v < u){
              if((Model_Params$B[u,v] == 0) & Results$B_estimates$B_final_pars$B[u,v] == 0){
                TNRB = TNRB + 1/(sum(Model_Params$B == 0) - meta_param$s*(meta_param$s-1)/2)
              }
              if((Model_Params$B[u,v] != 0) & Results$B_estimates$B_final_pars$B[u,v] != 0){
                TPRB = TPRB + 1/(sum(Model_Params$B != 0) - meta_param$s)
              }
            }
          }
        }
        for(u in (meta_param$s+1):meta_param$q){
          for(v in 1:meta_param$s){
            if((Model_Params$B[u,v] == 0) & Results$B_estimates$B_final_pars$B[u,v] == 0){
              TNRB = TNRB + 1/(sum(Model_Params$B == 0) - meta_param$s*(meta_param$s-1)/2)
            }
            if((Model_Params$B[u,v] != 0) & Results$B_estimates$B_final_pars$B[u,v] != 0){
              TPRB = TPRB + 1/(sum(Model_Params$B != 0) - meta_param$s)
            }
          }
        }
        BMSElist <- (2/(meta_param$s*(2*meta_param$q-meta_param$s-1)))*(Results$B_estimates$B_final_pars$B - Model_Params$B)^2
        itb1[which(i == c(50, 60, 70, 80))] <- itb1[which(i == c(50, 60, 70, 80))] + 1
      }
      if(optlatA == meta_param$m){
        TNRA <- 0
        TPRA <- 0
        for(u in 1:meta_param$m){
          for(v in 1:meta_param$m){
            if(v < u){
              if((Model_Params$A[u,v] == 0) & Results$A_estimates$A_final_pars$A[u,v] == 0){
                TNRA = TNRA + 1/(sum(Model_Params$A == 0) - meta_param$m*(meta_param$m-1)/2)
              }
              if((Model_Params$A[u,v] != 0) & Results$A_estimates$A_final_pars$A[u,v] != 0){
                TPRA = TPRA + 1/(sum(Model_Params$A != 0) - meta_param$m)
              }
            }
          }
        }
        for(u in (meta_param$m+1):meta_param$p){
          for(v in 1:meta_param$m){
            if((Model_Params$A[u,v] == 0) & Results$A_estimates$A_final_pars$A[u,v] == 0){
              TNRA = TNRA + 1/(sum(Model_Params$A == 0) - meta_param$m*(meta_param$m-1)/2)
            }
            if((Model_Params$A[u,v] != 0) & Results$A_estimates$A_final_pars$A[u,v] != 0){
              TPRA = TPRA + 1/(sum(Model_Params$A != 0) - meta_param$m)
            }
          }
        }
        AMSElist <- (2/(meta_param$m*(2*meta_param$q-meta_param$m-1)))*(Results$A_estimates$A_final_pars$A - Model_Params$A)^2
        ita1[which(i == c(50, 60, 70, 80))] <- ita1[which(i == c(50, 60, 70, 80))] + 1
      }
      if((optlatB == meta_param$s) & (optlatA == meta_param$m)){
        GMSElist <- (1/(meta_param$s*meta_param$m))*(Results$A_estimates$A_final_pars$Gamma - Model_Params$Gamma)^2
        itg1[which(i == c(50, 60, 70, 80))] <- itg1[which(i == c(50, 60, 70, 80))] + 1
      }
      iter <- iter + 1
      ResTable[[1]][iter,] <- c(optlatB,i,datfiles[k])
      ResTable[[2]][iter,] <- c(optlatA,i,datfiles[k])
      ResTable[[3]][iter,] <- c(opttunB,i,datfiles[k])
      ResTable[[4]][iter,] <- c(opttunA,i,datfiles[k])
      ResTable[[5]][iter,] <- c(sum(BMSElist),i,datfiles[k])
      ResTable[[6]][iter,] <- c(sum(AMSElist),i,datfiles[k])
      ResTable[[7]][iter,] <- c(sum(GMSElist),i,datfiles[k])
      ResTable[[8]][iter,] <- c(TNRB,i,datfiles[k])
      ResTable[[9]][iter,] <- c(TPRB,i,datfiles[k])
      ResTable[[10]][iter,] <- c(TNRA,i,datfiles[k])
      ResTable[[11]][iter,] <- c(TPRA,i,datfiles[k])
    }
  }
}

for(i in 1:11){
  colnames(ResTable[[i]]) <- c("Outcome", "n", "repnames")
  ResTable[[i]] <- transform(ResTable[[i]], Outcome = as.numeric(Outcome), n = as.numeric(n))
  ResTable[[i]] <- na.omit(ResTable[[i]])
  ResTable[[i]] <- ResTable[[i]][ResTable[[i]]$repnames != "",]
}



dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==50,-3], k = 5)
abline(h = 0.008, lty = 2)
db50 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==50,-3], eps = 0.008, MinPts = 5)
bout50 <- ResTable[[5]][ResTable[[5]]$n==50,][!db50$isseed,]

dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==60,-3], k = 5)
abline(h = 0.003, lty = 2)
db60 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==60,-3], eps = 0.003, MinPts = 5)
bout60 <- ResTable[[5]][ResTable[[5]]$n==60,][!db60$isseed,]

dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==70,-3], k = 5)
abline(h = 0.002, lty = 2)
db70 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==70,-3], eps = 0.002, MinPts = 5)
bout70 <- ResTable[[5]][ResTable[[5]]$n==70,][!db70$isseed,]

dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==80,-3], k = 5)
abline(h = 0.002, lty = 2)
db80 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==80,-3], eps = 0.002, MinPts = 5)
bout80 <- ResTable[[5]][ResTable[[5]]$n==80,][!db80$isseed,]



dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==50,-3], k = 5)
abline(h = 15, lty = 2)
db50 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==50,-3], eps = 15, MinPts = 5)
aout50 <- ResTable[[6]][ResTable[[6]]$n==50,][!db50$isseed,]

dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==60,-3], k = 5)
abline(h = 12, lty = 2)
db60 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==60,-3], eps = 12, MinPts = 5)
aout60 <- ResTable[[6]][ResTable[[6]]$n==60,][!db60$isseed,]

dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==70,-3], k = 5)
abline(h = 2, lty = 2)
db70 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==70,-3], eps = 2, MinPts = 5)
aout70 <- ResTable[[6]][ResTable[[6]]$n==70,][!db70$isseed,]

dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==80,-3], k = 5)
abline(h = 1, lty = 2)
db80 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==80,-3], eps = 1, MinPts = 5)
aout80 <- ResTable[[6]][ResTable[[6]]$n==80,][!db80$isseed,]



dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==50,-3], k = 5)
abline(h = 0.4, lty = 2)
db50 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==50,-3], eps = 0.4, MinPts = 5)
gout50 <- ResTable[[7]][ResTable[[7]]$n==50,][!db50$isseed,]

dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==60,-3], k = 5)
abline(h = 0.4, lty = 2)
db60 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==60,-3], eps = 0.4, MinPts = 5)
gout60 <- ResTable[[7]][ResTable[[7]]$n==60,][!db60$isseed,]

dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==70,-3], k = 5)
abline(h = 0.3, lty = 2)
db70 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==70,-3], eps = 0.3, MinPts = 5)
gout70 <- ResTable[[7]][ResTable[[7]]$n==70,][!db70$isseed,]

dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==80,-3], k = 5)
abline(h = 0.3, lty = 2)
db80 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==80,-3], eps = 0.3, MinPts = 5)
gout80 <- ResTable[[7]][ResTable[[7]]$n==80,][!db80$isseed,]



boutliers <- c(bout50$repnames, bout60$repnames, bout70$repnames, bout80$repnames)
aoutliers <- c(aout50$repnames, aout60$repnames, aout70$repnames, aout80$repnames)
goutliers <- c(gout50$repnames, gout60$repnames, gout70$repnames, gout80$repnames)

outliers <- replicate(11,c(),simplify=F)
outliers[[5]] <- boutliers
outliers[[6]] <- aoutliers
outliers[[7]] <- goutliers

AvgTables <- replicate(11,matrix(0,4,2),simplify=F)
for(i in 1:11){
  ResTable[[i]] <- ResTable[[i]][!ResTable[[i]][,3] %in% outliers[[i]],]
  AvgTables[[i]] <- aggregate(ResTable[[i]]$Outcome, list(ResTable[[i]]$n), FUN = mean, na.rm=TRUE, na.action=NULL)
  AvgTables[[i]] <- AvgTables[[i]][AvgTables[[i]][,1] != 0,]
}

propb <- itb1/totreps
propa <- ita1/totreps
propg <- itg1/totreps