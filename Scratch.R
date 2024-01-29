rm(list = ls())
library(mgcv)
library(dplyr)
library(lubridate)
library(viridis)
library(ggplot2)
library(tidyr)
library(ggcorrplot)


# Load data, drift GPS, and noise levels pre-process, and wind lease area

# Set the directory where your CSV files are located
csv_directory <- "F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023"

# Get a list of CSV files in the directory (adjust the pattern if needed)
csv_files <- list.files(path = csv_directory, 
                        pattern = "*.csv", full.names = TRUE)

# Combine the list of dataframes into a single dataframe
GPSdf <- do.call(rbind, lapply(csv_files, read.csv))

# Add UTC coords with the sp package
cord.dec = SpatialPoints(cbind(GPSdf$Longitude,GPSdf$Latitude),
                         proj4string=CRS("+proj=longlat"))
cord.dec = spTransform(cord.dec, CRS("+init=epsg:32610"))
GPSdf$UTMx =  coordinates(cord.dec)[,1]
GPSdf$UTMy =  coordinates(cord.dec)[,2]

# Convert to datetime, create date hour column
GPSdf$UTC=as.POSIXct(GPSdf$UTC,tz = 'UTC')
GPSdf$UTCDatehour = GPSdf$UTC
minute(GPSdf$UTCDatehour)<-0
second(GPSdf$UTCDatehour)<-0

# 
# # Add noise level data
# csv_directory='F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023 Noise Files'
# csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
# dataframes_list <- list()
# 
# # Loop through each CSV file load, change name
# for (csv_file in csv_files) {
#   # Read the CSV file
#   df <- read.csv(csv_file, header = TRUE)
#   
#   colnames(df)[1]<-'UTC'
#   
#   # Extract the first 10 characters from the filename
#   file_name <- substr(basename(csv_file), 1, 10)
#   
#   # Create a 'DriftName' column with the extracted filename
#   df$DriftName <- file_name
#   
#   # Append the dataframe to the list
#   dataframes_list[[file_name]] <- df
# }
# 
# 
# # Combine all dataframes into a single dataframe
# noiseDf <- bind_rows(dataframes_list)
# save(noiseDf, file = 'AllNL')

load('AllNL')


# Clean out data for drifts that don't have gps or noise levels
GPSdf= subset(GPSdf, DriftName %in% noiseDf$DriftName)
noiseDf= subset(noiseDf, DriftName %in% GPSdf$DriftName)


noise500 = noiseDf[, c(1,8, 25)]
noise20k = noiseDf[, c(1,24, 25)]

noise500reshape <- noise500 %>%
  spread(key = DriftName, value = TOL_500)


noise20kreshape <- noise20k %>%
  spread(key = DriftName, value = TOL_20000)

corrMat20khz <- cor(noise20kreshape[,c(2:8)], use = "pairwise.complete.obs") 
corrMat500hz <- cor(noise500reshape[,c(2:8)], use = "pairwise.complete.obs") 
pmat20 =cor_pmat(corrMat20khz)
# handling missing values if any


ggcorrplot(corrMat500hz, method = "circle", 
           lab = TRUE,  # to display correlation coefficients
           type = "upper",  # display lower half of matrix
           colors = c("#6D9EC1", "white", "red3"),
           title = 'Correlation 500 Hz')  # color scale

ggcorrplot(corrMat20khz, method = "circle", 
           lab = TRUE,  # to display correlation coefficients
           type = "upper",  # display lower half of matrix
           colors = c("#6D9EC1", "white", "red3"),
           title = 'Correlation 20 kHz')  # color scale

colnames(noise20k)[2]<-"NL"
colnames(noise500)[2]<-"NL"
noise20k$Band = '20kHz'
noise500$Band = '500Hz'

allNl = rbind(noise500, noise20k)
allNl$Band <- factor(allNl$Band , 
                     levels = c("500Hz", "20kHz"))


ggplot(allNl)+
  geom_line(aes(x = UTC, y=NL, group = DriftName, color = DriftName))+
  facet_wrap(~Band,nrow = 2)

