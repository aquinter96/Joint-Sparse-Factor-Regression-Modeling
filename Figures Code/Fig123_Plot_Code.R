library(lavaan)
library(lessSEM)
library(regsem)
library(ggplot2)

ResTable <- replicate(4,data.frame("n" = double(), "Matrix" = character(), "MSE" = as.double(), "Model" = character(), "Subgroup" = character()),simplify=F)

set.seed(2435)

source(file = "Model1_Meta_Param.R")
source(file = "Model1_Matrix_Param.R")

q <- meta_param$q
p <- meta_param$p
s <- meta_param$s
m <- meta_param$m

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("Lavrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])[[1]]
    ResTable[[1]][counter, 1] <- i
    ResTable[[1]][counter, 2] <- "B"
    ResTable[[1]][counter, 3] <- 1/(q*s-s*(s+1)/2)*sum((Model_Params$B - lavInspect(Output, what = "est")$lambda[1:q, 1:s])^2)
    ResTable[[1]][counter, 4] <- "Lavaan"
    ResTable[[1]][counter, 5] <- "Blavaan"
    ResTable[[1]][counter+1, 1] <- i
    ResTable[[1]][counter+1, 2] <- "A"
    ResTable[[1]][counter+1, 3] <- 1/(p*m-m*(m+1)/2)*sum((Model_Params$A - lavInspect(Output, what = "est")$lambda[(q+(1:p)), (s+(1:m))])^2)
    ResTable[[1]][counter+1, 4] <- "Lavaan"
    ResTable[[1]][counter+1, 5] <- "Alavaan"
    ResTable[[1]][counter+2, 1] <- i
    ResTable[[1]][counter+2, 2] <- "Gamma"
    ResTable[[1]][counter+2, 3] <- 1/(s*m)*sum((Model_Params$Gamma - lavInspect(Output, what = "est")$beta[(s+(1:m)), 1:s])^2)
    ResTable[[1]][counter+2, 4] <- "Lavaan"
    ResTable[[1]][counter+2, 5] <- "Gammalavaan" 
    counter <- counter + 3
  }
}

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("Lessrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])[[1]]
    Output <- lessSEM2Lavaan(Output, criterion = "BIC")
    ResTable[[2]][counter, 1] <- i
    ResTable[[2]][counter, 2] <- "B"
    ResTable[[2]][counter, 3] <- 1/(q*s-s*(s+1)/2)*sum((Model_Params$B - lavInspect(Output, what = "est")$lambda[1:q, 1:s])^2)
    ResTable[[2]][counter, 4] <- "LessSEM"
    ResTable[[2]][counter, 5] <- "Blesssem"
    ResTable[[2]][counter+1, 1] <- i
    ResTable[[2]][counter+1, 2] <- "A"
    ResTable[[2]][counter+1, 3] <- 1/(p*m-m*(m+1)/2)*sum((Model_Params$A - lavInspect(Output, what = "est")$lambda[(q+(1:p)), (s+(1:m))])^2)
    ResTable[[2]][counter+1, 4] <- "LessSEM"
    ResTable[[2]][counter+1, 5] <- "Alesssem"
    ResTable[[2]][counter+2, 1] <- i
    ResTable[[2]][counter+2, 2] <- "Gamma"
    ResTable[[2]][counter+2, 3] <- 1/(s*m)*sum((Model_Params$Gamma - lavInspect(Output, what = "est")$beta[(s+(1:m)), 1:s])^2)
    ResTable[[2]][counter+2, 4] <- "LessSEM"
    ResTable[[2]][counter+2, 5] <- "Gammalesssem"
    counter <- counter + 3
  }
}

regB <- matrix(0, nrow = q, ncol = s)
diag(regB[1:s, 1:s]) <- 1
regA <- matrix(0, nrow = p, ncol = m)
diag(regB[1:m, 1:m]) <- 1
regG <- matrix(0, nrow = m, ncol = s)

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("regrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])[[1]]
    
    regB[2:q,1] <- Output$final_pars[1:19]
    regB[3:q,2] <- Output$final_pars[20:37]
    regB[4:q,3] <- Output$final_pars[38:54]
    
    regA[2:p,1] <- Output$final_pars[55:73]
    regA[3:p,2] <- Output$final_pars[74:91]
    
    regG[1,] <- Output$final_pars[92:94]
    regG[2,] <- Output$final_pars[95:97]
    
    ResTable[[3]][counter, 1] <- i
    ResTable[[3]][counter, 2] <- "B"
    ResTable[[3]][counter, 3] <- 1/(q*s-s*(s+1)/2)*sum((Model_Params$B - regB)^2)
    ResTable[[3]][counter, 4] <- "Regsem"
    ResTable[[3]][counter, 5] <- "Bregsem"
    ResTable[[3]][counter+1, 1] <- i
    ResTable[[3]][counter+1, 2] <- "A"
    ResTable[[3]][counter+1, 3] <- 1/(p*m-m*(m+1)/2)*sum((Model_Params$A - regA)^2)
    ResTable[[3]][counter+1, 4] <- "Regsem"
    ResTable[[3]][counter+1, 5] <- "Aregsem"
    ResTable[[3]][counter+2, 1] <- i
    ResTable[[3]][counter+2, 2] <- "Gamma"
    ResTable[[3]][counter+2, 3] <- 1/(s*m)*sum((Model_Params$Gamma - regG)^2)
    ResTable[[3]][counter+2, 4] <- "Regsem"
    ResTable[[3]][counter+2, 5] <- "Gammaregsem"
    counter <- counter + 3
  }
}

counter <- 1
for(i in c(200, 400, 600, 800)){
  datfiles <- list.files(pattern=paste0("JSFRMrep.*n",i,".rds"))
  for(j in 1:length(datfiles)){
    Output <-readRDS(datfiles[j])
    ResTable[[4]][counter, 1] <- i
    ResTable[[4]][counter, 2] <- "B"
    ResTable[[4]][counter, 3] <- 1/(q*s-s*(s+1)/2)*sum((Model_Params$B - Output$B_estimates$B_final_pars$B)^2)
    ResTable[[4]][counter, 4] <- "JSFRM"
    ResTable[[4]][counter, 5] <- "Bjsfrm"
    ResTable[[4]][counter+1, 1] <- i
    ResTable[[4]][counter+1, 2] <- "A"
    ResTable[[4]][counter+1, 3] <- 1/(p*m-m*(m+1)/2)*sum((Model_Params$A - Output$A_estimates$A_final_pars$A)^2)
    ResTable[[4]][counter+1, 4] <- "JSFRM"
    ResTable[[4]][counter+1, 5] <- "Ajsfrm"
    ResTable[[4]][counter+2, 1] <- i
    ResTable[[4]][counter+2, 2] <- "Gamma"
    ResTable[[4]][counter+2, 3] <- 1/(s*m)*sum((Model_Params$Gamma - Output$A_estimates$A_final_pars$Gamma)^2)
    ResTable[[4]][counter+2, 4] <- "JSFRM"
    ResTable[[4]][counter+2, 5] <- "Gammajsfrm"
    counter <- counter + 3
  }
}

ResTable <- do.call(rbind, ResTable)
ResTable <- ResTable[complete.cases(ResTable),]

BMSEs <- ResTable[ResTable$Matrix == "B",]
AMSEs <- ResTable[ResTable$Matrix == "A",]
GMSEs <- ResTable[ResTable$Matrix == "Gamma",]

BMSEavg <- aggregate(BMSEs$MSE, list(BMSEs$n,BMSEs$Model), mean)
names(BMSEavg) <- c("n", "Model", "MSE")
AMSEavg <- aggregate(AMSEs$MSE, list(AMSEs$n,AMSEs$Model), mean)
names(AMSEavg) <- c("n", "Model", "MSE")
GMSEavg <- aggregate(GMSEs$MSE, list(GMSEs$n,GMSEs$Model), mean)
names(GMSEavg) <- c("n", "Model", "MSE")

ggplot(data = BMSEavg, aes(x = as.factor(n), y = log10(MSE), group = Model, color = Model)) +  geom_point() + geom_line() + theme_bw() + scale_fill_grey() + xlab("Sample Size") + ylab("log(MSE)")
ggplot(data = AMSEavg, aes(x = as.factor(n), y = log10(MSE), group = Model, color = Model)) +  geom_point() + geom_line() + theme_bw() + scale_fill_grey() + xlab("Sample Size") + ylab("log(MSE)")
ggplot(data = GMSEavg, aes(x = as.factor(n), y = log10(MSE), group = Model, color = Model)) +  geom_point() + geom_line() + theme_bw() + scale_fill_grey() + xlab("Sample Size") + ylab("log(MSE)")

BMSEs$rowname <- rownames(BMSEs)
AMSEs$rowname <- rownames(AMSEs)
GMSEs$rowname <- rownames(GMSEs)

Blav <- BMSEs[BMSEs$Model=="Lavaan",-c(2, 4, 5)]
dbscan::kNNdistplot(Blav[(Blav$n==200),-3], k = 5)
abline(h = 0.01, lty = 2)
db200 <- fpc::dbscan(Blav[Blav$n==200,-3], eps = 0.01, MinPts = 5)
if(!is.null(db200$isseed)){
  boutlav200 <- Blav[Blav$n==200,][!db200$isseed,]
}else{
  boutlav200 <- NULL
}

dbscan::kNNdistplot(Blav[(Blav$n==400),-3], k = 5)
abline(h = 0.005, lty = 2)
db400 <- fpc::dbscan(Blav[Blav$n==400,-3], eps = 0.005, MinPts = 5)
if(!is.null(db400$isseed)){
  boutlav400 <- Blav[Blav$n==400,][!db400$isseed,]
}else{
  boutlav400 <- NULL
}

dbscan::kNNdistplot(Blav[(Blav$n==600),-3], k = 5)
abline(h = 0.005, lty = 2)
db600 <- fpc::dbscan(Blav[Blav$n==600,-3], eps = 0.005, MinPts = 5)
if(!is.null(db600$isseed)){
  boutlav600 <- Blav[Blav$n==600,][!db600$isseed,]
}else{
  boutlav600 <- NULL
}

dbscan::kNNdistplot(Blav[(Blav$n==800),-3], k = 5)
abline(h = 0.005, lty = 2)
db800 <- fpc::dbscan(Blav[Blav$n==800,-3], eps = 0.005, MinPts = 5)
if(!is.null(db800$isseed)){
  boutlav800 <- Blav[Blav$n==800,][!db800$isseed,]
}else{
  boutlav800 <- NULL
}


Bless <- BMSEs[BMSEs$Model=="LessSEM",-c(2, 4, 5)]
dbscan::kNNdistplot(Bless[(Bless$n==200),-3], k = 5)
abline(h = 0.01, lty = 2)
db200 <- fpc::dbscan(Bless[Bless$n==200,-3], eps = 0.01, MinPts = 5)
if(!is.null(db200$isseed)){
  boutless200 <- Bless[Bless$n==200,][!db200$isseed,]
}else{
  boutless200 <- NULL
}

dbscan::kNNdistplot(Bless[(Bless$n==400),-3], k = 5)
abline(h = 0.01, lty = 2)
db400 <- fpc::dbscan(Bless[Bless$n==400,-3], eps = 0.01, MinPts = 5)
if(!is.null(db400$isseed)){
  boutless400 <- Bless[Bless$n==400,][!db400$isseed,]
}else{
  boutless400 <- NULL
}

dbscan::kNNdistplot(Bless[(Bless$n==600),-3], k = 5)
abline(h = 0.005, lty = 2)
db600 <- fpc::dbscan(Bless[Bless$n==600,-3], eps = 0.005, MinPts = 5)
if(!is.null(db600$isseed)){
  boutless600 <- Bless[Bless$n==600,][!db600$isseed,]
}else{
  boutless600 <- NULL
}

dbscan::kNNdistplot(Bless[(Bless$n==800),-3], k = 5)
abline(h = 0.005, lty = 2)
db800 <- fpc::dbscan(Bless[Bless$n==800,-3], eps = 0.005, MinPts = 5)
if(!is.null(db800$isseed)){
  boutless800 <- Bless[Bless$n==800,][!db800$isseed,]
}else{
  boutless800 <- NULL
}


Breg <- BMSEs[BMSEs$Model=="Regsem",-c(2, 4, 5, 6)]
dbscan::kNNdistplot(Breg[(Breg$n==200),-3], k = 5)
abline(h = 0.01, lty = 2)
db200 <- fpc::dbscan(Breg[Breg$n==200,-3], eps = 0.01, MinPts = 5)
if(!is.null(db200$isseed)){
  boutreg200 <- Breg[Breg$n==200,][!db200$isseed,]
}else{
  boutreg200 <- NULL
}

dbscan::kNNdistplot(Breg[(Breg$n==400),-3], k = 5)
abline(h = 0.001, lty = 2)
db400 <- fpc::dbscan(Breg[Breg$n==400,-3], eps = 0.001, MinPts = 5)
if(!is.null(db400$isseed)){
  boutreg400 <- Breg[Breg$n==400,][!db400$isseed,]
}else{
  boutreg400 <- NULL
}

dbscan::kNNdistplot(Breg[(Breg$n==600),-3], k = 5)
abline(h = 0.0003, lty = 2)
db600 <- fpc::dbscan(Breg[Breg$n==600,-3], eps = 0.0003, MinPts = 5)
if(!is.null(db600$isseed)){
  boutreg600 <- Breg[Breg$n==600,][!db600$isseed,]
}else{
  boutreg600 <- NULL
}

dbscan::kNNdistplot(Breg[(Breg$n==800),-3], k = 5)
abline(h = 0.0005, lty = 2)
db800 <- fpc::dbscan(Breg[Breg$n==800,-3], eps = 0.0005, MinPts = 5)
if(!is.null(db800$isseed)){
  boutreg800 <- Breg[Breg$n==800,][!db800$isseed,]
}else{
  boutreg800 <- NULL
}


Bjsf <- BMSEs[BMSEs$Model=="JSFRM",-c(2, 4, 5)]
dbscan::kNNdistplot(Bjsf[(Bjsf$n==200),-3], k = 5)
abline(h = 0.005, lty = 2)
db200 <- fpc::dbscan(Bjsf[Bjsf$n==200,-3], eps = 0.005, MinPts = 5)
if(!is.null(db200$isseed)){
  boutjsf200 <- Bjsf[Bjsf$n==200,][!db200$isseed,]
}else{
  boutjsf200 <- NULL
}

dbscan::kNNdistplot(Bjsf[(Bjsf$n==400),-3], k = 5)
abline(h = 0.001, lty = 2)
db400 <- fpc::dbscan(Bjsf[Bjsf$n==400,-3], eps = 0.001, MinPts = 5)
if(!is.null(db400$isseed)){
  boutjsf400 <- Bjsf[Bjsf$n==400,][!db400$isseed,]
}else{
  boutjsf400 <- NULL
}

dbscan::kNNdistplot(Bjsf[(Bjsf$n==600),-3], k = 5)
abline(h = 0.0008, lty = 2)
db600 <- fpc::dbscan(Bjsf[Bjsf$n==600,-3], eps = 0.0008, MinPts = 5)
if(!is.null(db600$isseed)){
  boutjsf600 <- Bjsf[Bjsf$n==600,][!db600$isseed,]
}else{
  boutjsf600 <- NULL
}

dbscan::kNNdistplot(Bjsf[(Bjsf$n==800),-3], k = 5)
abline(h = 0.0003, lty = 2)
db800 <- fpc::dbscan(Bjsf[Bjsf$n==800,-3], eps = 0.0003, MinPts = 5)
if(!is.null(db800$isseed)){
  boutjsf800 <- Bjsf[Bjsf$n==800,][!db800$isseed,]
}else{
  boutjsf800 <- NULL
}




Alav <- AMSEs[AMSEs$Model=="Lavaan",-c(2, 4, 5)]
dbscan::kNNdistplot(Alav[(Alav$n==200),-3], k = 5)
abline(h = 0.1*10^8, lty = 2)
db200 <- fpc::dbscan(Alav[Alav$n==200,-3], eps = 0.1*10^8, MinPts = 5)
if(!is.null(db200$isseed)){
  aoutlav200 <- Alav[Alav$n==200,][!db200$isseed,]
}else{
  aoutlav200 <- NULL
}

dbscan::kNNdistplot(Alav[(Alav$n==400),-3], k = 5)
abline(h = 1.3*10^7, lty = 2)
db400 <- fpc::dbscan(Alav[Alav$n==400,-3], eps = 1.3*10^7, MinPts = 5)
if(!is.null(db400$isseed)){
  aoutlav400 <- Alav[Alav$n==400,][!db400$isseed,]
}else{
  aoutlav400 <- NULL
}

dbscan::kNNdistplot(Alav[(Alav$n==600),-3], k = 5)
abline(h = 10^7, lty = 2)
db600 <- fpc::dbscan(Alav[Alav$n==600,-3], eps = 10^7, MinPts = 5)
if(!is.null(db600$isseed)){
  aoutlav600 <- Alav[Alav$n==600,][!db600$isseed,]
}else{
  aoutlav600 <- NULL
}

dbscan::kNNdistplot(Alav[(Alav$n==800),-3], k = 5)
abline(h = 10^5, lty = 2)
db800 <- fpc::dbscan(Alav[Alav$n==800,-3], eps = 10^5, MinPts = 5)
if(!is.null(db800$isseed)){
  aoutlav800 <- Alav[Alav$n==800,][!db800$isseed,]
}else{
  aoutlav800 <- NULL
}


Aless <- AMSEs[AMSEs$Model=="LessSEM",-c(2, 4, 5)]
dbscan::kNNdistplot(Aless[(Aless$n==200),-3], k = 5)
abline(h = 8*10^7, lty = 2)
db200 <- fpc::dbscan(Aless[Aless$n==200,-3], eps = 8*10^7, MinPts = 5)
if(!is.null(db200$isseed)){
  aoutless200 <- Aless[Aless$n==200,][!db200$isseed,]
}else{
  aoutless200 <- NULL
}

dbscan::kNNdistplot(Aless[(Aless$n==400),-3], k = 5)
abline(h = 1.1*10^7, lty = 2)
db400 <- fpc::dbscan(Aless[Aless$n==400,-3], eps = 1.1*10^7, MinPts = 5)
if(!is.null(db400$isseed)){
  aoutless400 <- Aless[Aless$n==400,][!db400$isseed,]
}else{
  aoutless400 <- NULL
}

dbscan::kNNdistplot(Aless[(Aless$n==600),-3], k = 5)
abline(h = 10^7, lty = 2)
db600 <- fpc::dbscan(Aless[Aless$n==600,-3], eps = 10^7, MinPts = 5)
if(!is.null(db600$isseed)){
  aoutless600 <- Aless[Aless$n==600,][!db600$isseed,]
}else{
  aoutless600 <- NULL
}

dbscan::kNNdistplot(Aless[(Aless$n==800),-3], k = 5)
abline(h = 1.5*10^5, lty = 2)
db800 <- fpc::dbscan(Aless[Aless$n==800,-3], eps = 1.5*10^5, MinPts = 5)
if(!is.null(db800$isseed)){
  aoutless800 <- Aless[Aless$n==800,][!db800$isseed,]
}else{
  aoutless800 <- NULL
}


Areg <- AMSEs[AMSEs$Model=="Regsem",-c(2, 4, 5)]
dbscan::kNNdistplot(Areg[(Areg$n==200),-3], k = 5)
abline(h = 0.001, lty = 2)
db200 <- fpc::dbscan(Areg[Areg$n==200,-3], eps = 0.001, MinPts = 5)
if(!is.null(db200$isseed)){
  aoutreg200 <- Areg[Areg$n==200,][!db200$isseed,]
}else{
  aoutreg200 <- NULL
}

dbscan::kNNdistplot(Areg[(Areg$n==400),-3], k = 5)
abline(h = 0.001, lty = 2)
db400 <- fpc::dbscan(Areg[Areg$n==400,-3], eps = 0.001, MinPts = 5)
if(!is.null(db400$isseed)){
  aoutreg400 <- Areg[Areg$n==400,][!db400$isseed,]
}else{
  aoutreg400 <- NULL
}

dbscan::kNNdistplot(Areg[(Areg$n==600),-3], k = 5)
abline(h = 0.0006, lty = 2)
db600 <- fpc::dbscan(Areg[Areg$n==600,-3], eps = 0.0006, MinPts = 5)
if(!is.null(db600$isseed)){
  aoutreg600 <- Areg[Areg$n==600,][!db600$isseed,]
}else{
  aoutreg600 <- NULL
}

dbscan::kNNdistplot(Areg[(Areg$n==800),-3], k = 5)
abline(h = 0.0002, lty = 2)
db800 <- fpc::dbscan(Areg[Areg$n==800,-3], eps = 0.0002, MinPts = 5)
if(!is.null(db800$isseed)){
  aoutreg800 <- Areg[Areg$n==800,][!db800$isseed,]
}else{
  aoutreg800 <- NULL
}


Ajsf <- AMSEs[AMSEs$Model=="JSFRM",-c(2, 4, 5)]
dbscan::kNNdistplot(Ajsf[(Ajsf$n==200),-3], k = 5)
abline(h = 1, lty = 2)
db200 <- fpc::dbscan(Ajsf[Ajsf$n==200,-3], eps = 1, MinPts = 5)
if(!is.null(db200$isseed)){
  aoutjsf200 <- Ajsf[Ajsf$n==200,][!db200$isseed,]
}else{
  aoutjsf200 <- NULL
}

dbscan::kNNdistplot(Ajsf[(Ajsf$n==400),-3], k = 5)
abline(h = 1, lty = 2)
db400 <- fpc::dbscan(Ajsf[Ajsf$n==400,-3], eps = 1, MinPts = 5)
if(!is.null(db400$isseed)){
  aoutjsf400 <- Ajsf[Ajsf$n==400,][!db400$isseed,]
}else{
  aoutjsf400 <- NULL
}

dbscan::kNNdistplot(Ajsf[(Ajsf$n==600),-3], k = 5)
abline(h = 0.0002, lty = 2)
db600 <- fpc::dbscan(Ajsf[Ajsf$n==600,-3], eps = 0.0002, MinPts = 5)
if(!is.null(db600$isseed)){
  aoutjsf600 <- Ajsf[Ajsf$n==600,][!db600$isseed,]
}else{
  aoutjsf600 <- NULL
}

dbscan::kNNdistplot(Ajsf[(Ajsf$n==800),-3], k = 5)
abline(h = 0.0002, lty = 2)
db800 <- fpc::dbscan(Ajsf[Ajsf$n==800,-3], eps = 0.0002, MinPts = 5)
if(!is.null(db800$isseed)){
  aoutjsf800 <- Ajsf[Ajsf$n==800,][!db800$isseed,]
}else{
  aoutjsf800 <- NULL
}




Glav <- GMSEs[GMSEs$Model=="Lavaan",-c(2, 4, 5)]
dbscan::kNNdistplot(Glav[(Glav$n==200),-3], k = 5)
abline(h = 100, lty = 2)
db200 <- fpc::dbscan(Glav[Glav$n==200,-3], eps = 100, MinPts = 5)
if(!is.null(db200$isseed)){
  goutlav200 <- Glav[Glav$n==200,][!db200$isseed,]
}else{
  goutlav200 <- NULL
}

dbscan::kNNdistplot(Glav[(Glav$n==400),-3], k = 5)
abline(h = 0.01, lty = 2)
db400 <- fpc::dbscan(Glav[Glav$n==400,-3], eps = 0.01, MinPts = 5)
if(!is.null(db400$isseed)){
  goutlav400 <- Glav[Glav$n==400,][!db400$isseed,]
}else{
  goutlav400 <- NULL
}

dbscan::kNNdistplot(Glav[(Glav$n==600),-3], k = 5)
abline(h = 0.001, lty = 2)
db600 <- fpc::dbscan(Glav[Glav$n==600,-3], eps = 0.001, MinPts = 5)
if(!is.null(db600$isseed)){
  goutlav600 <- Glav[Glav$n==600,][!db600$isseed,]
}else{
  goutlav600 <- NULL
}

dbscan::kNNdistplot(Glav[(Glav$n==800),-3], k = 5)
abline(h = 0.01, lty = 2)
db800 <- fpc::dbscan(Glav[Glav$n==800,-3], eps = 0.01, MinPts = 5)
if(!is.null(db800$isseed)){
  goutlav800 <- Glav[Glav$n==800,][!db800$isseed,]
}else{
  goutlav800 <- NULL
}


Gless <- GMSEs[GMSEs$Model=="LessSEM",-c(2, 4, 5)]
dbscan::kNNdistplot(Gless[(Gless$n==200),-3], k = 5)
abline(h = 100, lty = 2)
db200 <- fpc::dbscan(Gless[Gless$n==200,-3], eps = 100, MinPts = 5)
if(!is.null(db200$isseed)){
  goutless200 <- Gless[Gless$n==200,][!db200$isseed,]
}else{
  goutless200 <- NULL
}

dbscan::kNNdistplot(Gless[(Gless$n==400),-3], k = 5)
abline(h = 0.01, lty = 2)
db400 <- fpc::dbscan(Gless[Gless$n==400,-3], eps = 0.01, MinPts = 5)
if(!is.null(db400$isseed)){
  goutless400 <- Gless[Gless$n==400,][!db400$isseed,]
}else{
  goutless400 <- NULL
}

dbscan::kNNdistplot(Gless[(Gless$n==600),-3], k = 5)
abline(h = 0.005, lty = 2)
db600 <- fpc::dbscan(Gless[Gless$n==600,-3], eps = 0.005, MinPts = 5)
if(!is.null(db600$isseed)){
  goutless600 <- Gless[Gless$n==600,][!db600$isseed,]
}else{
  goutless600 <- NULL
}

dbscan::kNNdistplot(Gless[(Gless$n==800),-3], k = 5)
abline(h = 0.01, lty = 2)
db800 <- fpc::dbscan(Gless[Gless$n==800,-3], eps = 0.01, MinPts = 5)
if(!is.null(db800$isseed)){
  goutless800 <- Gless[Gless$n==800,][!db800$isseed,]
}else{
  goutless800 <- NULL
}


Greg <- GMSEs[GMSEs$Model=="Regsem",-c(2, 4, 5)]
dbscan::kNNdistplot(Greg[(Greg$n==200),-3], k = 5)
abline(h = 0.01, lty = 2)
db200 <- fpc::dbscan(Greg[Greg$n==200,-3], eps = 0.01, MinPts = 5)
if(!is.null(db200$isseed)){
  goutreg200 <- Greg[Greg$n==200,][!db200$isseed,]
}else{
  goutreg200 <- NULL
}

dbscan::kNNdistplot(Greg[(Greg$n==400),-3], k = 5)
abline(h = 0.003, lty = 2)
db400 <- fpc::dbscan(Greg[Greg$n==400,-3], eps = 0.003, MinPts = 5)
if(!is.null(db400$isseed)){
  goutreg400 <- Greg[Greg$n==400,][!db400$isseed,]
}else{
  goutreg400 <- NULL
}

dbscan::kNNdistplot(Greg[(Greg$n==600),-3], k = 5)
abline(h = 0.0006, lty = 2)
db600 <- fpc::dbscan(Greg[Greg$n==600,-3], eps = 0.0006, MinPts = 5)
if(!is.null(db600$isseed)){
  goutreg600 <- Greg[Greg$n==600,][!db600$isseed,]
}else{
  goutreg600 <- NULL
}

dbscan::kNNdistplot(Greg[(Greg$n==800),-3], k = 5)
abline(h = 0.001, lty = 2)
db800 <- fpc::dbscan(Greg[Greg$n==800,-3], eps = 0.001, MinPts = 5)
if(!is.null(db800$isseed)){
  goutreg800 <- Greg[Greg$n==800,][!db800$isseed,]
}else{
  goutreg800 <- NULL
}


Gjsf <- GMSEs[GMSEs$Model=="JSFRM",-c(2, 4, 5)]
dbscan::kNNdistplot(Gjsf[(Gjsf$n==200),-3], k = 5)
abline(h = 0.01, lty = 2)
db200 <- fpc::dbscan(Gjsf[Gjsf$n==200,-3], eps = 0.01, MinPts = 5)
if(!is.null(db200$isseed)){
  goutjsf200 <- Gjsf[Gjsf$n==200,][!db200$isseed,]
}else{
  goutjsf200 <- NULL
}

dbscan::kNNdistplot(Gjsf[(Gjsf$n==400),-3], k = 5)
abline(h = 0.01, lty = 2)
db400 <- fpc::dbscan(Gjsf[Gjsf$n==400,-3], eps = 0.01, MinPts = 5)
if(!is.null(db400$isseed)){
  goutjsf400 <- Gjsf[Gjsf$n==400,][!db400$isseed,]
}else{
  goutjsf400 <- NULL
}

dbscan::kNNdistplot(Gjsf[(Gjsf$n==600),-3], k = 5)
abline(h = 0.0005, lty = 2)
db600 <- fpc::dbscan(Gjsf[Gjsf$n==600,-3], eps = 0.0005, MinPts = 5)
if(!is.null(db600$isseed)){
  goutjsf600 <- Gjsf[Gjsf$n==600,][!db600$isseed,]
}else{
  goutjsf600 <- NULL
}

dbscan::kNNdistplot(Gjsf[(Gjsf$n==800),-3], k = 5)
abline(h = 0.0004, lty = 2)
db800 <- fpc::dbscan(Gjsf[Gjsf$n==800,-3], eps = 0.0004, MinPts = 5)
if(!is.null(db800$isseed)){
  goutjsf800 <- Gjsf[Gjsf$n==800,][!db800$isseed,]
}else{
  goutjsf800 <- NULL
}


boutliers <- c(boutlav200$rowname,boutlav400$rowname,boutlav600$rowname,boutlav800$rowname,
               boutless200$rowname,boutless400$rowname,boutless600$rowname,boutless800$rowname,
               boutreg200$rowname,boutreg400$rowname,boutreg600$rowname,boutreg800$rowname,
               boutjsf200$rowname,boutjsf400$rowname,boutjsf600$rowname,boutjsf800$rowname)
aoutliers <- c(aoutlav200$rowname,aoutlav400$rowname,aoutlav600$rowname,aoutlav800$rowname,
               aoutless200$rowname,aoutless400$rowname,aoutless600$rowname,aoutless800$rowname,
               aoutreg200$rowname,aoutreg400$rowname,aoutreg600$rowname,aoutreg800$rowname,
               aoutjsf200$rowname,aoutjsf400$rowname,aoutjsf600$rowname,aoutjsf800$rowname)
goutliers <- c(goutlav200$rowname,goutlav400$rowname,goutlav600$rowname,goutlav800$rowname,
               goutless200$rowname,goutless400$rowname,goutless600$rowname,goutless800$rowname,
               goutreg200$rowname,goutreg400$rowname,goutreg600$rowname,goutreg800$rowname,
               goutjsf200$rowname,goutjsf400$rowname,goutjsf600$rowname,goutjsf800$rowname)

BMSEnew <- BMSEs[!BMSEs$rowname %in% boutliers,]
AMSEnew <- AMSEs[!AMSEs$rowname %in% aoutliers,]
GMSEnew <- GMSEs[!GMSEs$rowname %in% goutliers,]

BMSEnewavg <- aggregate(BMSEnew$MSE, list(BMSEnew$n,BMSEnew$Model), mean)
names(BMSEnewavg) <- c("n", "Model", "MSE")
AMSEnewavg <- aggregate(AMSEnew$MSE, list(AMSEnew$n,AMSEnew$Model), mean)
names(AMSEnewavg) <- c("n", "Model", "MSE")
GMSEnewavg <- aggregate(GMSEnew$MSE, list(GMSEnew$n,GMSEnew$Model), mean)
names(GMSEnewavg) <- c("n", "Model", "MSE")

ggplot(data = BMSEnewavg, aes(x = as.factor(n), y = log10(MSE), group = Model, color = Model)) +  geom_point() + geom_line() + theme_bw() + scale_fill_grey() + xlab("Sample Size") + ylab("log(MSE)")
ggplot(data = AMSEnewavg, aes(x = as.factor(n), y = log10(MSE), group = Model, color = Model)) +  geom_point() + geom_line() + theme_bw() + scale_fill_grey() + xlab("Sample Size") + ylab("log(MSE)")
ggplot(data = GMSEnewavg, aes(x = as.factor(n), y = log10(MSE), group = Model, color = Model)) +  geom_point() + geom_line() + theme_bw() + scale_fill_grey() + xlab("Sample Size") + ylab("log(MSE)")
