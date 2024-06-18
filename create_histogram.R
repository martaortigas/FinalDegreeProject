# R Script to build the histogram displaying the frequency of values of dxy_WNL_DSE – dxy_WSL_DSE from Table S1.

library(readxl)

rm(list=ls())

setwd("D:/Internship_TFG/documentos_prácticas")

dataset <- read_excel("Dxy_WSL_introgr_to_DSE_AIob.xlsx")

View(dataset)

x <- hist(dataset$`dxy_WNL_DSE - dxy_WSL_DSE`, main = "dxy_WNL_DSE - dxy_WSL_DSE", xlab = "dxy_WNL_DSE-dxy_WSL_DSE")

mean_value <- mean(dataset$`dxy_WNL_DSE - dxy_WSL_DSE`)
max_value <- max(dataset$`dxy_WNL_DSE - dxy_WSL_DSE`)
min_value <- min(dataset$`dxy_WNL_DSE - dxy_WSL_DSE`)

abline(v=mean_value, col="red", lwd=2, lty=2)  
abline(v=max_value, col="green", lwd=2, lty=2) 
abline(v=min_value, col="blue", lwd=2, lty=2)  

# Add a legend
legend_x <- mean(x$breaks) * 0.1  # You can adjust this as needed
legend_y <- max(x$counts) * 1  # You can adjust this as needed

legend(legend_x, legend_y, 
       legend=c(paste("Mean: ", round(mean_value, 2)),
                paste("Max: ", round(max_value, 2)),
                paste("Min: ", round(min_value, 5))),
       col=c("red", "green", "blue"), lty=2, lwd=2)
