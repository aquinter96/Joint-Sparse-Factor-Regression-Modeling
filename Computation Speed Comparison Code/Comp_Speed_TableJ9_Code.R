
JSFRMtime <- matrix(NA, nrow = 200, ncol = 4)
lavtime <- matrix(NA, nrow = 200, ncol = 4)
regtime <- matrix(NA, nrow = 200, ncol = 4)
lesstime <- matrix(NA, nrow = 200, ncol = 4)

iter <- 1
for(j in c(200, 400, 600, 800)){
  for(k in c("JSFRM", "Lav", "Less", "reg")){
    iter <- 1
    datfiles <- list.files(pattern=paste0(k, "rep.*n", j, ".rds"))
    for(l in 1:length(datfiles)){
      if(file.exists(datfiles[l])){
        Output <- readRDS(datfiles[l])
        if(k == "Lav"){
          lavtime[iter, which(c(200, 400, 600, 800) == j)] <- Output[[2]][1]
        }
        if(k == "reg"){
          regtime[iter, which(c(200, 400, 600, 800) == j)] <- Output[[2]][1]
        }
        if(k == "Less"){
          lesstime[iter, which(c(200, 400, 600, 800) == j)] <- Output[[2]][1]
        }
        if(k == "JSFRM"){
          JSFRMtime[iter, which(c(200, 400, 600, 800) == j)] <- Output$tot.time[1]
        }
        iter <- iter + 1
      }
    }
  }
}
