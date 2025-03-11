library(lavaan)
library(gdata)
library(readr)

args <- commandArgs(TRUE)

Country <- as.character(args[1])

Country <- gsub("_", " ", Country)

Country.Results <- readRDS(paste0(Country, ".Results.rds"))

s_hat <- Country.Results$B_estimates$s_opt
m_hat <- Country.Results$A_estimates$m_opt

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

Xnewnames <- rep("", 27*s_hat)

for(i in 1:s_hat){
  for(j in 1:27){
    Xnewnames[j+(i - 1)*27] <- paste0("x", j+(i - 1)*27, "*", colnames(X.Country.train)[j])
  }
}

Xvals <- rep("", s_hat)

for(i in 1:s_hat){
  Xvals[i] <- paste(c(Xnewnames[(1+(i - 1)*27):(i*27)]), collapse = " + ")
}

for(i in 1:s_hat) {
  assign(paste0("X", i, "model"), paste(paste0(sprintf("W%d =~ ", i), Xvals[i]), sep = '\n'))
}

Xmodel <- c()

Xstrings <- c()

for(i in 1:s_hat) {
  mv(from = paste0("X", i, "model"), "Xstrings")
  Xmodel <- c(Xmodel, Xstrings)
}

Wvals <- substr(paste(sprintf("W%d +", 1:s_hat), collapse = ' '), 1 ,nchar(paste(sprintf("W%d +", 1:s_hat), collapse = ' '))  - 2)
Zmodel <- paste(paste0(sprintf("Z%d ~ ", 1:m_hat), Wvals), sep = '\n')
Wvar <- paste(sprintf("W%d ~~ 1*W%d", 1:s_hat, 1:s_hat), sep = '\n')
Wcov <- c()

k <- 1

Xfree <- c()
for(i in 1:(27*s_hat)){
  Xfree[i] <- paste0("x", i)
}

Xfixed <- c()

for(i in 1:s_hat){
  Xmodel[i] <- gsub(paste0("x", i+27*(i-1), "[*]", colnames(X.Country.train)[i], " "), paste0("1*", colnames(X.Country.train)[i], " "), Xmodel[i])
  Xfixed[i] <- i+27*(i-1)
  for(j in 1:s_hat){
    if(i < j){
      Wcov[k] <- paste(sprintf("W%d ~~ 0*W%d", i, j), sep = '\n')
      k <- k+1
    }
  }
}

if(s_hat > 1){
  for(i in 2:s_hat){
    for(j in 1:(i-1)){
      Xmodel[i] <- gsub(paste0("x", j+27*(i - 1)), paste0("0"), Xmodel[i])
      Xfixed <- c(Xfixed, j+27*(i - 1))
    }
  }
}

Xfree <- Xfree[-Xfixed]

Ynewnames <- rep("", 22*m_hat)

for(i in 1:m_hat){
  for(j in 1:22){
    Ynewnames[j+(i - 1)*22] <- paste0("y", j+(i - 1)*22, "*", colnames(Y.Country.train)[j])
  }
}

Yvals <- rep("", m_hat)

for(i in 1:m_hat){
  Yvals[i] <- paste(c(Ynewnames[(1+(i - 1)*22):(i*22)]), collapse = " + ")
}

for(i in 1:m_hat) {
  assign(paste0("Y", i, "model"), paste(paste0(sprintf("Z%d =~ ", i), Yvals[i]), sep = '\n'))
}

Ymodel <- c()

Ystrings <- c()

for(i in 1:m_hat) {
  mv(from = paste0("Y", i, "model"), "Ystrings")
  Ymodel <- c(Ymodel, Ystrings)
}

k <- 1

Yfree <- c()
for(i in 1:(22*m_hat)){
  Yfree[i] <- paste0("y", i)
}

Yfixed <- c()

for(i in 1:m_hat){
  Ymodel[i] <- gsub(paste0("y", i+22*(i-1), "[*]", colnames(Y.Country.train)[i], " "), paste0("1*", colnames(Y.Country.train)[i], " "), Ymodel[i])
  Yfixed[i] <- i+22*(i-1)
}

if(m_hat > 1){
  for(i in 2:m_hat){
    for(j in 1:(i-1)){
      Ymodel[i] <- gsub(paste0("y", j+22*(i - 1)), paste0("0"), Ymodel[i])
      Yfixed <- c(Yfixed, j+22*(i - 1))
    }
  }
}

Yfree <- Yfree[-Yfixed]

set.seed(2435)

lavmodel <- c(Xmodel, Ymodel, Zmodel, Wvar, Wcov)
lavdat <- cbind(X.Country.train, Y.Country.train)

lav_mod <- sem(lavmodel, data = lavdat, meanstructure = T, orthogonal = T)

saveRDS(lav_mod, file=paste(Country,".lavaan.Results.rds", sep=""))