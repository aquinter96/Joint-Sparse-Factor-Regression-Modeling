library(ggplot2)

setwd("~/Paper1_SIM_Real")

Countries <- c("Australia", "Austria", "Brazil", "Canada", "Croatia", "Germany", "Greece",
               "Indonesia", "Ireland", "Malaysia", "Netherlands", "Pakistan", "Panama", "Philippines",
               "Poland", "Portugal", "Serbia", "Slovakia", "South Korea", "Spain",
               "Switzerland", "Turkey", "United Kingdom", "United States")

Country.plot <- NULL

for(Country in Countries){
  
  test = readRDS(paste0(Country, ".Results.rds"))
  mydat <- data.frame(1:20, test$B_estimates$BICs)
  names(mydat) <- c("s_hat", "XBICs")
  
  Country.plot <- ggplot(mydat, aes(x=as.factor(s_hat), y=XBICs, fill=as.factor(s_hat))) + geom_bar(stat = "identity") + theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +  geom_text(aes(label=round(XBICs)), position=position_dodge(width=0.9), vjust=-0.25) + scale_x_discrete(name ="Value of s") + ggtitle(paste0(Country, " X Model BICs"))
  
  print(Country.plot)
  
  ggsave(filename = paste0(Country, "_XBICs.png"),Country.plot, path = "~/Paper1_SIM_Real")
  
}

Country.plot <- NULL

for(Country in Countries){
  
  test = readRDS(paste0(Country, ".Alt_Results.rds"))
  mydat <- data.frame(1:20, test$A_estimates$BICs)
  names(mydat) <- c("m_hat", "YBICs")
  
  print(paste0(Country, ":", test$A_estimates$m_opt))
  
  Country.plot <- ggplot(mydat, aes(x=as.factor(m_hat), y=YBICs, fill=as.factor(m_hat))) + geom_bar(stat = "identity") + theme_minimal() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +  geom_text(aes(label=round(YBICs)), position=position_dodge(width=0.9), vjust=-0.25) + scale_x_discrete(name ="Value of m") + ggtitle(paste0(Country, " Y Model BICs"))
  
  print(Country.plot)
  
  ggsave(filename = paste0(Country, "_YBICs.png"),Country.plot, path = "~/Paper1_SIM_Real")
  
}