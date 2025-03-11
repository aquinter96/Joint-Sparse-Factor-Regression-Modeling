#Model1_Table_Code.R

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

source(file = "Model1_Meta_Param.R")
source(file = "Model1_Matrix_Param.R")

totreps <- rep(0, 4)

for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("Model1rep.*n",i,".rds"))
  totreps[which(i == c(200, 400, 600, 800))] <- length(datfiles)
  print(length(datfiles))
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
        itb1[which(i == c(200, 400, 600, 800))] <- itb1[which(i == c(200, 400, 600, 800))] + 1
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
        ita1[which(i == c(200, 400, 600, 800))] <- ita1[which(i == c(200, 400, 600, 800))] + 1
      }
      if((optlatB == meta_param$s) & (optlatA == meta_param$m)){
        GMSElist <- (1/(meta_param$s*meta_param$m))*(Results$A_estimates$A_final_pars$Gamma - Model_Params$Gamma)^2
        itg1[which(i == c(200, 400, 600, 800))] <- itg1[which(i == c(200, 400, 600, 800))] + 1
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
      if(!is.na(sum(BMSElist)) & sum(BMSElist) > 10){
        print(c("B", sum(BMSElist), datfiles[k]))
        itb2[1] <- itb2[1] + 1
      }
      if(!is.na(sum(AMSElist)) & sum(AMSElist) > 10){
        print(c("A", sum(AMSElist), datfiles[k]))
        ita2[1] <- ita2[1] + 1
      }
      if(!is.na(sum(GMSElist)) & sum(GMSElist) > 10){
        print(c("G", sum(GMSElist), datfiles[k]))
        itg2[1] <- itg2[1] + 1
      }
    }
  }
}

for(i in 1:11){
  colnames(ResTable[[i]]) <- c("Outcome", "n", "repnames")
  ResTable[[i]] <- transform(ResTable[[i]], Outcome = as.numeric(Outcome), n = as.numeric(n))
  ResTable[[i]] <- na.omit(ResTable[[i]])
  ResTable[[i]] <- ResTable[[i]][ResTable[[i]]$repnames != "",]
}



dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==200,-3], k = 5)
abline(h = 0.006, lty = 2)
db200 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==200,-3], eps = 0.006, MinPts = 5)
bout200 <- ResTable[[5]][ResTable[[5]]$n==200,][!db200$isseed,]

dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==400,-3], k = 5)
abline(h = 0.001, lty = 2)
db400 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==400,-3], eps = 0.001, MinPts = 5)
bout400 <- ResTable[[5]][ResTable[[5]]$n==400,][!db400$isseed,]

dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==600,-3], k = 5)
abline(h = 0.001, lty = 2)
db600 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==600,-3], eps = 0.001, MinPts = 5)
bout600 <- ResTable[[5]][ResTable[[5]]$n==600,][!db600$isseed,]

dbscan::kNNdistplot(ResTable[[5]][ResTable[[5]]$n==800,-3], k = 5)
abline(h = 0.0005, lty = 2)
db800 <- fpc::dbscan(ResTable[[5]][ResTable[[5]]$n==800,-3], eps = 0.0005, MinPts = 5)
bout800 <- ResTable[[5]][ResTable[[5]]$n==800,][!db800$isseed,]



dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==200,-3], k = 5)
abline(h = 0.005, lty = 2)
db200 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==200,-3], eps = 0.005, MinPts = 5)
aout200 <- ResTable[[6]][ResTable[[6]]$n==200,][!db200$isseed,]

dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==400,-3], k = 5)
abline(h = 0.0003, lty = 2)
db400 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==400,-3], eps = 0.0003, MinPts = 5)
aout400 <- ResTable[[6]][ResTable[[6]]$n==400,][!db400$isseed,]

dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==600,-3], k = 5)
abline(h = 0.0002, lty = 2)
db600 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==600,-3], eps = 0.0002, MinPts = 5)
aout600 <- ResTable[[6]][ResTable[[6]]$n==600,][!db600$isseed,]

dbscan::kNNdistplot(ResTable[[6]][ResTable[[6]]$n==800,-3], k = 5)
abline(h = 0.0001, lty = 2)
db800 <- fpc::dbscan(ResTable[[6]][ResTable[[6]]$n==800,-3], eps = 0.0001, MinPts = 5)
aout800 <- ResTable[[6]][ResTable[[6]]$n==800,][!db800$isseed,]



dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==200,-3], k = 5)
abline(h = 0.008, lty = 2)
db200 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==200,-3], eps = 0.008, MinPts = 5)
gout200 <- ResTable[[7]][ResTable[[7]]$n==200,][!db200$isseed,]

dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==400,-3], k = 5)
abline(h = 0.0006, lty = 2)
db400 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==400,-3], eps = 0.0006, MinPts = 5)
gout400 <- ResTable[[7]][ResTable[[7]]$n==400,][!db400$isseed,]

dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==600,-3], k = 5)
abline(h = 0.0005, lty = 2)
db600 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==600,-3], eps = 0.0005, MinPts = 5)
gout600 <- ResTable[[7]][ResTable[[7]]$n==600,][!db600$isseed,]

dbscan::kNNdistplot(ResTable[[7]][ResTable[[7]]$n==800,-3], k = 5)
abline(h = 0.0004, lty = 2)
db800 <- fpc::dbscan(ResTable[[7]][ResTable[[7]]$n==800,-3], eps = 0.0004, MinPts = 5)
gout800 <- ResTable[[7]][ResTable[[7]]$n==800,][!db800$isseed,]



boutliers <- c(bout200$repnames, bout400$repnames, bout600$repnames, bout800$repnames)
aoutliers <- c(aout200$repnames, aout400$repnames, aout600$repnames, aout800$repnames)
goutliers <- c(gout200$repnames, gout400$repnames, gout600$repnames, gout800$repnames)

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