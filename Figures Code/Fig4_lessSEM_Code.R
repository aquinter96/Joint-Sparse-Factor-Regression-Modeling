library(lavaan)
library(stringr)
library(lessSEM)

args <- commandArgs(TRUE)
datset <- as.numeric(args[1]) 
sampsize <- as.numeric(args[2])

set.seed(2435)

source(file = "Model1_Meta_Param.R")
source(file = "Model1_Matrix_Param.R")

source(file = "Data_Generation.R")

set.seed(datset)

dat <- Data_Generation(n = sampsize, Model_Params=Model_Params)

Split.ind <- sample(seq_len(nrow(dat$X)), size = floor(0.60*nrow(dat$X)))
X.train <- dat$X[Split.ind,]
X.test <- dat$X[-Split.ind,]
Y.train <- dat$Y[Split.ind,]
Y.test <- dat$Y[-Split.ind,]

colnames(X.train) <- sprintf("X%d", 1:q)
colnames(Y.train) <- sprintf("Y%d", 1:p)

mydat <- cbind(X.train, Y.train)

X1vals <- substr(paste(sprintf("x%d*X%d +", 1:q, 1:q), collapse = ' '), 1 ,nchar(paste(sprintf("x%d*X%d +", 1:q, 1:q), collapse = ' '))  - 2)
X2vals <- substr(paste(sprintf("x%d*X%d +", (q+1):(2*q), 1:q), collapse = ' '), 1 ,nchar(paste(sprintf("x%d*X%d +", (q+1):(2*q), 1:q), collapse = ' '))  - 2)
X3vals <- substr(paste(sprintf("x%d*X%d +", (2*q+1):(3*q), 1:q), collapse = ' '), 1 ,nchar(paste(sprintf("x%d*X%d +", (2*q+1):(3*q), 1:q), collapse = ' '))  - 2)
Y1vals <- substr(paste(sprintf("y%d*Y%d +", 1:p, 1:p), collapse = ' '), 1 ,nchar(paste(sprintf("y%d*Y%d +", 1:p, 1:p), collapse = ' '))  - 2)
Y2vals <- substr(paste(sprintf("y%d*Y%d +", (p+1):(2*p), 1:p), collapse = ' '), 1 ,nchar(paste(sprintf("y%d*Y%d +", (p+1):(2*p), 1:p), collapse = ' '))  - 2)
Wvals <- substr(paste(sprintf("w%d*W%d +", 1:s, 1:s), collapse = ' '), 1 ,nchar(paste(sprintf("w%d*W%d +", 1:s, 1:s), collapse = ' '))  - 2)
X1model <- paste(paste0(sprintf("W%d =~ ", 1), X1vals), sep = '\n')
X2model <- paste(paste0(sprintf("W%d =~ ", 2), X2vals), sep = '\n')
X3model <- paste(paste0(sprintf("W%d =~ ", 3), X3vals), sep = '\n')
Y1model <- paste(paste0(sprintf("Z%d =~ ", 1), Y1vals), sep = '\n')
Y2model <- paste(paste0(sprintf("Z%d =~ ", 2), Y2vals), sep = '\n')
Xmodel <- c(X1model, X2model, X3model)
Ymodel <- c(Y1model, Y2model)
Wvar <- paste(sprintf("W%d ~~ 1*W%d", 1:s, 1:s), sep = '\n')
Wcov <- c()
Zvar <- paste(sprintf("Z%d ~~ 1*Z%d", 1:m, 1:m), sep = '\n')
Zcov <- c()

k <- 1

for(i in 1:s){
  Xmodel[i] <- gsub(paste0("x", i+q*(i-1), "[*]X", i, " "), paste0("1*X", i, " "), Xmodel[i])
  for(j in 1:s){
    if(i < j){
      Wcov[k] <- paste(sprintf("W%d ~~ 0*W%d", i, j), sep = '\n')
      k <- k+1
    }
  }
}
Xmodel[2] <- gsub(paste0("x21"), paste0("0"), Xmodel[2])
Xmodel[3] <- gsub(paste0("x41"), paste0("0"), Xmodel[3])
Xmodel[3] <- gsub(paste0("x42"), paste0("0"), Xmodel[3])
k <- 1
for(i in 1:m){
  Ymodel[i] <- gsub(paste0("y", i+p*(i-1), "[*]Y", i, " "), paste0("1*Y", i, " "), Ymodel[i])
  for(j in 1:m){
    if(i < j){
      Zcov[k] <- paste(sprintf("Z%d ~~ 0*Z%d", i, j), sep = '\n')
      k <- k+1
    }
  }
}
Ymodel[2] <- gsub(paste0("y21"), paste0("0"), Ymodel[2])
Zmodel <- c("Z1 ~ W1 + W2 + W3", "Z2 ~ W1 + W2 + W3")

mymodel <- c(Xmodel, Ymodel, Zmodel, Wvar)
lav_mod <- sem(mymodel, data = mydat, meanstructure = T, orthogonal = T)

free_params <- parTable(lav_mod)$label
free_params <- free_params[nzchar(free_params)]

less_mod <- adaptiveLasso(lav_mod, regularized = free_params, lambdas = seq(0, 15, 0.1))

saveRDS(less_mod, file=paste0("Lessrep", datset, "n", sampsize, ".rds"))