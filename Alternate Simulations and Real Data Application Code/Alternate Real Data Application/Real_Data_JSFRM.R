library(readr)
library(parallel)

args <- commandArgs(TRUE)
Country <- as.character(args[1])

Country <- gsub("_", " ", Country)

COVIDiSTRESS_global_survey_May_30_2020_final_cleaned_file_ <- read_csv("COVIDiSTRESS global survey May 30 2020 (___final cleaned file___).csv")

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
Country.train <- NULL
Country.train$X <- as.matrix(X.Country.train)
Country.train$Y <- as.matrix(Y.Country.train)

source(file = "Lasso.R")
source(file = "B_inits.R")
source(file = "Estep_X.R")
source(file = "Mstep_X.R")
source(file = "logLik_X.R")
source(file = "Convergence_check.R")
source(file = "BIC_X.R")
source(file = "EMAlgBAdLassoCV.R")
source(file = "OverallBAlg.R")

Singular_ErrorX <- function(tuningpB, Data, initial_Model_Params, weights){
  print(weights)
  return(tryCatch(EMAlgBAdLassoCV(tuningpB, Data, initial_Model_Params, weights), error = function(e) NULL))
}

numCores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
cl <- makeCluster(as.numeric(numCores))
clusterExport(cl, c("lasso", "Estep_X", "Mstep_X", "logLik_X", "Convergence_check", "BIC_X", "EMAlgBAdLassoCV"))

Country.Results.B <- OverallBAlg(Country.train$X, s_seq = 4)

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

Country.Results <- OverallAGAlg(Country.train, Country.Results.B, m_seq = 4)

stopCluster(cl)

saveRDS(Country.Results, file=paste0(Country,".Alt2_Results.rds"))