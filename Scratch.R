rm(list = ls())
library(mgcv)
library(dplyr)
library(lubridate)
library(viridis)
library(ggplot2)
library(tidyr)
library(ggcorrplot)
library(sp)
library(geosphere)


########################################################################
# Load data, drift GPS, and noise levels pre-process, and wind lease area
########################################################################
# GPS Directory
csv_directory <- "F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023"

# List of the csv files
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



# Add noise level data
csv_directory='F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023 Noise Files'
csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
dataframes_list <- list()

# Loop through each CSV file load, change name
for (csv_file in csv_files) {
  # Read the CSV file
  df <- read.csv(csv_file, header = TRUE)
  
  # Extract the first and eighth columns
  colnames(df)[1]<-'UTC'
  
  # Extract the first 10 characters from the filename
  file_name <- substr(basename(csv_file), 1, 10)
  
  # Create a 'DriftName' column with the extracted filename
  df$DriftName <- file_name
  
  # Append the dataframe to the list
  dataframes_list[[file_name]] <- df
}

# Combine all dataframes into a single dataframe
noiseDf <- bind_rows(dataframes_list)
noiseDf$UTC <- as.POSIXct(noiseDf$UTC, 
                          format = "%Y-%m-%dT%H:%M:%OSZ",
                          tz='UTC')




# Clean out data for drifts that don't have gps or noise levels
GPSdf= subset(GPSdf, DriftName %in% noiseDf$DriftName)
noiseDf= subset(noiseDf, DriftName %in% GPSdf$DriftName)


noise500 = noiseDf[, c(1,8, 25)]
noise20k = noiseDf[, c(1,24, 25)]

noise500reshape <- noise500 %>%
  spread(key = DriftName, value = TOL_500)

noise20kreshape <- noise20k %>%
  spread(key = DriftName, value = TOL_20000)

colnames(noise20k)[2]<-"NL"
colnames(noise500)[2]<-"NL"
noise20k$Band = '20kHz'
noise500$Band = '500Hz'

allNl = rbind(noise500, noise20k)
allNl$Band <- factor(allNl$Band , 
                     levels = c("500Hz", "20kHz"))
allNl$UTC<- as.POSIXct(allNl$UTC, 
                       format = "%Y-%m-%dT%H:%M:%OSZ",
                       tz='UTC')

# Add lat and lon to the noise levels

# Link the gps and the noise levels
allNl$Lon = NaN
allNl$Lat = NaN

DriftNames = unique(GPSdf$DriftName)
GPSdf$NL<- NaN

for(drift in DriftNames){
  GPSsub = subset(GPSdf, DriftName == drift)
  NLsub =  subset(allNl, DriftName==drift)
  
  UTMflon <- approxfun(GPSsub$UTC, GPSsub$Longitude)
  UTMflat <- approxfun(GPSsub$UTC, GPSsub$Latitude)
  
  ##############################################################
  # Calculate the range from the whale to the GPS, TDOA and RL
  ############################################################
  
  # Lat/lon/ of the drift when the call was produced
  allNl$Lon[allNl$DriftName==drift] =UTMflon(NLsub$UTC)
  allNl$Lat[allNl$DriftName==drift] = UTMflat(NLsub$UTC)
  
}

allNl = allNl[!is.na(allNl$Lat),]


########################################################################
# Speed and distance travelled by the drifters
########################################################################

library(geosphere)

# Get the speed between gps points
# Calculate speed
calculate_speed <- function(df) {
  df <- df %>%
    arrange(DriftName, UTC) %>%
    group_by(DriftName) %>%
    mutate(
      Distance = distVincentySphere(cbind(Lon1 = lag(Longitude),
                                          Lat1 = lag(Latitude)), 
                                    cbind(Lon2 = Longitude, Lat2 = Latitude)),
      TimeDiff = as.numeric(difftime(UTC, lag(UTC), units = "secs")),
      Speed = Distance / TimeDiff
    ) %>%
    select(-TimeDiff)
  
  return(df)
}

GPSdf <- calculate_speed(GPSdf)


# Function to calculate total distance for each drifter
calculate_total_distance <- function(df) {
  summary_df <- df %>%
    arrange(DriftName, UTC) %>%
    group_by(DriftName) %>%
    summarise(
      TotalDistance = distVincentySphere(cbind(Lon1 = first(Longitude), Lat1 = first(Latitude)), 
                                         cbind(Lon2 = last(Longitude), Lat2 = last(Latitude)))
    )
  
  return(summary_df)
}

# Apply the function to your dataframe
summary_table <- calculate_total_distance(GPSdf)




########################################################################
# Explore temporal autocorrelation
########################################################################

corrMat20khz <- cor(noise20kreshape[,c(2:8)], use = "pairwise.complete.obs") 
corrMat500hz <- cor(noise500reshape[,c(2:8)], use = "pairwise.complete.obs") 


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



ggplot(allNl)+
  facet_wrap(~Band,nrow = 2)+
  geom_line(aes(x = UTC, y=NL, group = DriftName, color = DriftName))+
  scale_x_datetime(breaks = seq(from = min(allNl$UTC),
                                to = max(allNl$UTC), by = "14 hours"), 
                   labels = scales::date_format("%y-%m-%d")) +
  ylab(expression("Med. Thrid Ocave Noise Levels (dB"["rms"]*" )"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,  vjust = 0.5))



########################################################################
# Prep for Fields Analysis- calcualted median hourly levels
########################################################################

# Calculate the median hourly noise level for each DriftName, day, and hour
median_hourly_noise <- allNl %>%
  group_by(Band, Date = date(UTC), Hour = hour(UTC)) %>%
  summarise(medianNL = median(NL, na.rm = TRUE))

allNl$Date = date(allNl$UTC)
allNl$Hour = hour(allNl$UTC)

# Merge the median hourly noise back to the original data frame
allNl <- allNl %>%
  left_join(median_hourly_noise, by = c( "Date", "Hour", "Band"))

# Find differences in noise levels and replot
allNl$RelNoise = allNl$NL- allNl$medianNL

# Plot check
ggplot(allNl)+ 
  facet_wrap(~Band,nrow = 2)+
  geom_line(aes(x= UTC, y=RelNoise, color= DriftName))



# Look at the correlation matrix between normalized data
noise500 = noiseDf[, c(1,8, 25)]
noise20k = noiseDf[, c(1,24, 25)]

aa =allNl[allNl$Band=='500Hz', 
          c('UTC','RelNoise','DriftName')] 

normNoiseReshapeLF <- allNl[allNl$Band=='500Hz', 
                            c('UTC','RelNoise','DriftName')] %>%
  spread(key = DriftName, value = RelNoise)

normNoiseReshapeHF <- allNl[allNl$Band !='500Hz', 
                            c('UTC','RelNoise','DriftName')] %>%
  spread(key = DriftName, value = RelNoise)


corrMat20khz <- cor(normNoiseReshapeHF[,c(2:8)], use = "pairwise.complete.obs") 
corrMat500hz <- cor(normNoiseReshapeLF[,c(2:8)], use = "pairwise.complete.obs") 


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

# Ok, so by removing the median noise levels we have reduced the temporal
# autocorrelation significantly. We can see that there is an increasing 
# correlation with more closely spased instruments and hadrly any with 
# instruments further apart. This holds for both low and high frequency values


# Summarize by hour, because this is crazy big
# Calculate the median hourly noise level for each DriftName, day, and hour
simplifNL <- allNl %>%
  group_by(Band, DriftName, Date = date(UTC), Hour = hour(UTC)) %>%
  summarise(NL = median(NL, na.rm = TRUE),
            Lat = median(Lat, na.rm= TRUE),
            Lon = median(Lon, na.rm= TRUE),
            RelNoise = median(RelNoise, na.rm =TRUE),
            UTC = median(UTC, na.rm =TRUE))


#######################################################################
# Fields Package
#######################################################################


library(fields)

# Clear out any NA values

simplifNL <- simplifNL[complete.cases(simplifNL), , drop = FALSE]
loc<-cbind(simplifNL$Lon[simplifNL$Band== '500Hz'],
           simplifNL$Lat[simplifNL$Band== '500Hz']) #locations


tseries = simplifNL$UTC[simplifNL$Band== '500Hz'] 
telapsed = (tseries-min(tseries))

nl.lf = simplifNL$NL[simplifNL$Band== '500Hz']
nl.hf = simplifNL$NL[simplifNL$Band!= '500Hz']

nl.norm.lf = simplifNL$RelNoise[simplifNL$Band== '500Hz']
nl.norm.hf = simplifNL$RelNoise[simplifNL$Band!= '500Hz']

## Data exploration using quilt plot ###
quilt.plot(loc, nl.lf, nx=100, ny=100)
quilt.plot(loc, nl.hf, nx=100, ny=100)


## Data exploration using quilt plot ###
quilt.plot(loc, nl.norm.lf, nx=100, ny=100)
quilt.plot(loc, nl.norm.hf, nx=100, ny=100)

# Spatial process
grd <- list(x = seq(from = min(loc[,1]), to = max(loc[,1]), length.out = 200),
            y= seq(from = min(loc[,2]), to = max(loc[,2]), length.out = 400))

# Lf spatial processes normalized and not
obj.lf.norm<- spatialProcess( loc, nl.norm.lf,)
out.p.lf.norm<-predictSurface.Krig( obj.lf.norm, 
                                    grid.list=grd, extrap=TRUE)

obj.hf.norm<- spatialProcess( loc, nl.norm.hf,
                              cov.args = list(Covariance="Exponential"))
out.p.hf.norm<-predictSurface.Krig( obj.hf.norm, 
                                    grid.list=grd, extrap=TRUE)



# Find the standard errors for the normalized plots
out.p.lf<-predictSurfaceSE(obj.lf.norm, 
                           grid.list=grd, extrap=TRUE) 
out.p.hf<-predictSurfaceSE( obj.hf.norm, 
                            grid.list=grd, extrap=TRUE) 

image.plot( out.p.lf, col=larry.colors())
image.plot( out.p.hf, col=larry.colors())

# These plots of the standard errors give us a good idea of the spatial extent
# of valid prediction. For Low frequency sounds, we still have errors less than
# 1 db for most of the habitat covered. For normalized HF values standard error
# is at or above 2db across the majority of the habitat. This is approximately what
# we expect

imagePlot( out.p.lf, col=larry.colors(), main = "Standard errors for LF NL values")
points(loc[,1], loc[,2], pch=19, cex=.5)
title("Standard errors for LF NL values")

imagePlot( out.p.hf, col=larry.colors(), main = "Standard errors for HF NL values")
points(loc[,1], loc[,2], pch=19, cex=.5)



###########################################################################
# Create a variogram to look at spatial correlation
###########################################################################

# We are looking at these plots to understand the range at wich we can reasonably
# say a prediction can be made. In other words, we are looking for high correlation
# values or low variance values such that we can say a value measured at x1, y1
# can be used to predict a value at x2, y2. This ability will decrease with 
# range and we are looking for the point at which that happens. 
# My current thinking is that this evaluation should be done with the raw data
# then the range cutoff should be used for the normalized data. 


# Raw data first
v.hf <- vgram(loc=loc, y=nl.hf, N=30,lon.lat=TRUE) #use 30 bins
boxplotVGram(v.hf, breaks=v.hf$breaks,lon.lat=TRUE, ylab="sqrt(Variance)",
             xlab="distance", main = 'Emperical Variogram Raw HF')

plot(v.hf$stats["mean",]~v.hf$centers,
     main='Binned Variogram Raw HF',
     ylab="sqrt(Variance)", xlab="distance")

# In the raw HF correlation we can a major peak around 30km and a lesser peak
# around 5km. Again, I would expect that the correlation at ~30km correlates
# with the onset of one or both of the torms


v.lf <- vgram(loc=loc, y=nl.lf, N=30,lon.lat=TRUE) #use 30 bins
boxplotVGram(v.lf, breaks = v.lf$breaks,lon.lat=TRUE, 
             ylab="sqrt(Variance)", xlab="distance", 
             main = 'Emperical Variogram Raw LF')

plot(v.lf$stats["mean",]~v.lf$centers,
     main='Binned Variogram Raw LF',
     ylab="sqrt(Variance)", xlab="distance")
# In the raw LF correlation we can see two peaks at ~25 km and ~ 35km. I think
# this is likely the period when the storms rolled through. 



## Now look at the normalized variograms. Theoretically this should take care
# of the storm activity


v.hf.norm <- vgram(loc=loc, y=nl.norm.hf, N=20,lon.lat=TRUE) #use 30 bins
boxplotVGram(v.hf.norm, breaks=v.hf.norm$breaks,lon.lat=TRUE, ylab="sqrt(Variance)",
             xlab="distance", main = 'Emperical Variogram Normalized HF')
plot(v.hf.norm$stats["mean",]~v.hf.norm$centers,main='Binned Variogram HF Normalized',
     ylab="sqrt(Variance)", xlab="distance")
# This is quite intersting, We have the variation decreasing and therefore
# the correlation in values increasing as a function of range between the sensors


# Here we get two vary distinct peaks in variation as a function of range
v.lf.norm <- vgram(loc=loc, y=nl.norm.lf, N=30,lon.lat=TRUE) #use 30 bins
boxplotVGram(v.lf.norm, breaks=v.lf.norm$breaks,lon.lat=TRUE, 
             ylab="sqrt(Variance)", xlab="distance", 
             main = 'Emperical Variogram Normalized LF')

plot(v.lf.norm$stats["mean",]~v.lf.norm$centers,main='Binned Variogram LF Normalized',
     ylab="sqrt(Variance)", xlab="distance")



# Plot covariance functions vs variograms to estimate
plot(v.lf$stats["mean",1:11]~v.lf$centers[1:11],main='Binned Semivariogram')



##providing initial guesses and bounds
ls.exponential <- function(par){ # par = (rho,aRange,sigma^2)
  theoretical.vgram <- par[1] * (1 - exp(-v.lf$centers[1:17] / par[2])) + par[3]
  sum( (theoretical.vgram - v.lf$stats["mean",1:17])^2 ) }

##providing initial guesses and bounds
out.exponential <- optim(par=c(2,50,1),
                         fn=ls.exponential,
                         method="L-BFGS", 
                         lower = c(0,0,0), upper = c(80,100,80))

rho.exponential <- out.exponential$par[1]
aRange.exponential <- out.exponential$par[2]
sigma2.exponential <- out.exponential$par[3]

plot(v.lf$stats["mean",1:17]~v.lf$centers[1:17], main="Exponential vs. Matern")
lines(c(sigma2.exponential + 
          rho.exponential*(1-exp(-v.lf$centers/aRange.exponential)))~v.lf$centers,col=1)

# Mattern plot isn't working so we will skip it.

##providing initial guesses and bounds
#Matern plot 


###############################################################################
# Using spatila process to fit spatial models, rather than step by step covariance
# selection
######################################################################

fit.lf.norm<- spatialProcess(x = loc, 
                             y= nl.norm.lf, 
                             Distance = "rdist.earth")

surface(fit.lf.norm)
title("Kriging Predictions")

fit.lf.norm.exp<-spatialProcess(x = loc,
                                y= nl.norm.lf,Distance = "rdist.earth",
                                cov.args = list(Covariance ="Exponential"))

surface(fit.lf.norm.exp)

# plot the fitted model to evaluate fit
set.panel(1,2)
plot(fit.lf.norm)
plot(fit.lf.norm.exp)


x<-cbind(allNl$Lon[allNl$Band== '500Hz'],
           allNl$Lat[allNl$Band== '500Hz']) #locations
y<- allNl$RelNoise[allNl$Band== '500Hz']

# Clean duplicated data
y =y[!duplicated(x)]
x =x[!duplicated(x),]

# 
# # I've used the lambda and arange values from the spatial
# # process fit
# out<- mKrig(x,y,
#             cov.function="stationary.cov",
#             Covariance="Matern",Distance = 'rdist.earth',
#             aRange =1.408, lambda = 2)
# 
# out1<- mKrig(x,y,
#             cov.function="stationary.cov",
#             Covariance="Matern",Distance = 'rdist.earth',
#             aRange =5, lambda = .2)
# 
# out2<- mKrig(x,y,
#              cov.function="stationary.cov",
#              Covariance="Matern",Distance = 'rdist.earth',
#              aRange =1.4, lambda = .2)
# 
out3<- mKrig(x,y,
             Covariance="Matern",Distance = 'rdist.earth',
             aRange =1.408, lambda = 2)
# 
# 
# out4<- mKrig(x,y,
#              Covariance="Matern",Distance = 'rdist.earth',
#              aRange =2, lambda = .2)
# 
# 
# out5<- mKrig(x,y,
#              Covariance="Matern",
#              Distance = 'rdist.earth',
#              aRange = 10, lambda = .1)
# 
# out6<- mKrig(x,y,
#              Covariance="Matern",
#              Distance = 'rdist.earth',
#              aRange = 10, lambda = .01)
# 
# 
# out7<- mKrig(x,y,
#              Covariance="Matern",
#              Distance = 'rdist.earth',
#              aRange = 10, lambda = 2)
# 
# 
# out8<- mKrig(x,y,
#              Covariance="Matern",
#              Distance = 'rdist.earth',
#              aRange = 8, lambda = 1)

# 
# surface(out)
# surface(out1)
# surface(out2)
# surface(out3) # lowest p-values
# surface(out4)
# surface(out5)
# surface(out6)
# surface(out7) # Best

# #########################################################################
# run k-fold cross validation
# ########################################################################

allNl$Prediction =NA

for(drift in DriftNames){
  
  kfoldData = fieldModelResults(NLData = allNl,
                               band = '500Hz', 
                               predDrift =drift)
  
  # Clean duplicated data
  x = kfoldData[[1]]
  y = kfoldData[[2]]
  y =y[!duplicated(x)]
  x =x[!duplicated(x),]
  
  out<- mKrig(x,y,
               Covariance="Matern",Distance = 'rdist.earth',
               aRange =1.408, lambda = 2)
  
  allNl$Prediction[allNl$DriftName== drift]<- predict(out, kfoldData[[3]])
  
}

allNl$PredictionError = abs(allNl$Prediction - allNl$RelNoise)


# #########################################################################
# # See if we can fit with full data- no ran for like an hour
# ########################################################################
# 
# aa = allNl[allNl$Band == '500Hz',]
# aa$loc = paste0(aa$Lat, aa$Lon)
# aa$dup = duplicated(aa$loc)
# aa = aa[aa$dup == FALSE,]
# 
# 
# 
# loc<-cbind(aa$Lon[aa$Band== '500Hz'],
#            aa$Lat[aa$Band== '500Hz']) #locations
# 
# tseries = aa$UTC[aa$Band== '500Hz']
# telapsed = (tseries-min(tseries))
# 
# nl.lf = aa$NL[aa$Band== '500Hz']
# nl.norm.lf = aa$RelNoise[aa$Band== '500Hz']
# 
# keepIdx = which(!is.na(loc[,1]))
# 
# loc= loc[keepIdx,]
# nl.norm.lf= nl.norm.lf[keepIdx]
# 
# obj2<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
#                        cov.args = list(Covariance ="Exponential"))
# ########################################################################

# 
# # exponential vs Wendland covariance function
# tt = as.numeric(telapsed-mean(telapsed))/60
# 
# obj.norm <- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth")
# 
# 
# obj2<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
#                        cov.args = list(Covariance ="Exponential")) #nope
# 
# obj3<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
#                        cov.args = list(Covariance = "Wendland",
#                                        dimension = 1, k = 1))
# 
# obj4<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
#                        cov.args = list(Covariance = "Matern"))
# 
# rbind(obj.norm$summary, obj2$summary, obj3$summary, obj3$summary)
# # Plot 1 shows data vs. predicted values, and plot 2 shows predicted values vs.
# # residuals. Plot 3 shows the criteria to select the smoothing parameter λ = σ
# # 2/ρ
# set.panel(1,2)
# plot(obj.norm)
# plot(obj2)
# plot(obj3)
# plot(obj4)
# 
# 
# obj.raw.tt <- spatialProcess( x = loc, y = nl.lf, Z = tt,
#                               Distance = "rdist.earth")
# obj2.raw.tt<- spatialProcess( loc, nl.lf, Distance = "rdist.earth",
#                               cov.args = list(Covariance ="Exponential"))#nope
# 
# obj3.raw.tt<- spatialProcess( loc,  y = nl.lf, Z = tt,
#                               Distance = "rdist.earth",
#                               cov.args = list(Covariance = "Wendland",
#                                               dimension = 1, k = 1))
# 
# obj4.raw.tt<- spatialProcess( loc, y = nl.lf, Z = tt,
#                               Distance = "rdist.earth",
#                               cov.args = list(Covariance = "Matern"))
# 
# 
# rbind(obj.raw.tt$summary,obj2.raw.tt$summary,obj3.raw.tt$summary,
#       obj4.raw.tt$summary)
# 
# set.panel(1,2)
# plot(obj.raw.tt)
# plot(obj2.raw.tt)
# plot(obj3.raw.tt)
# plot(obj4.raw.tt)#best
# 
# 
# # Samething but normalized
# obj.norm.tt <- spatialProcess( x = loc, y = nl.norm.lf, Z = tt,
#                                Distance = "rdist.earth")
# obj2.norm.tt<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
#                                cov.args = list(Covariance ="Exponential"))
# 
# obj3.norm.tt<- spatialProcess( loc,  y = nl.norm.lf, Z = tt,
#                                Distance = "rdist.earth",
#                                cov.args = list(Covariance = "Wendland",
#                                                dimension = 1, k = 1))
# 
# obj4.norm.tt<- spatialProcess( loc, y = nl.norm.lf, Z = tt,
#                                Distance = "rdist.earth",
#                                cov.args = list(Covariance = "Matern"))
# 
# 
# rbind(obj.norm.tt$summary,
#       obj2.norm.tt$summary,
#       obj3.norm.tt$summary,
#       obj4.norm.tt$summary)
# 
# set.panel(1,1)
# plot(obj.norm.tt)
# plot(obj2.norm.tt)
# plot(obj3.norm.tt)
# plot(obj4.norm.tt)
# 
# 
# # thin plate spline model

predgrd <- list(x = seq(from = min(loc[,1]-.25), to = max(loc[,1]+.25), length.out = 200),
                y= seq(from = min(loc[,2]-.25), to = max(loc[,2]+.25), length.out = 400))

 obj.TPS.raw <- Tps( loc, nl.lf, Z= tt,lon.lat = TRUE,m = 3)
 out.p <-predictSurface(obj.TPS.raw, grid.list=grd, ZGrid= predgrd, extrap=TRUE)
 
library(sp)
library(RColorBrewer)
image.plot( out.p, col=brewer.pal(9,"RdGy"))

###########################################################################
# Plot the best model
###########################################################################


# 
# obj4.raw.tt<- spatialProcess( loc, y = nl.lf, Z = s(tt,20),
#                               Distance = "rdist.earth",
#                               cov.args = list(Covariance = "Matern"))

predgrd <- list(x = seq(from = min(loc[,1]-.25), to = max(loc[,1]+.25), length.out = 200),
                y= seq(from = min(loc[,2]-.25), to = max(loc[,2]+.25), length.out = 400))



out.preds<-predictSurface.Krig(object = out7, grid.list = grd,
                               extrap=TRUE) 

image.plot( out.preds, col=larry.colors())
out.SEs<-predictSurfaceSE(object = out3, grid.list = grd, extrap=TRUE) 

# 
# out.preds<-predictSurface.Krig(object = obj4.raw.tt, grid.list = grd,
#                                ZGrid = predgrd, 
#                                extrap=TRUE) 
# 
# image.plot( out.preds, col=larry.colors())
# 
# out.SEs<-predictSurfaceSE(object = obj4.raw.tt, grid.list = grd,
#                           ZGrid = predgrd, 
#                           extrap=TRUE) 




###########################################################################
# Idiot check this with an animation
###########################################################################


library(gganimate)
library(ggplot2)


# Make multiple projections- minutes
tsteps = quantile(tt, probs = seq(0,1,by=.05))


dataOut = array()
x = array()
y= array()
t = array()

for(ii in 1:length(tsteps)){
  predgrd <- list(x = seq(from = min(loc[,1]), to = max(loc[,1]), length.out = 200),
                  y= seq(from = min(loc[,2]), to = max(loc[,2]), length.out = 400),
                  z = matrix(nrow = 200, ncol = 400, data = tsteps[ii]))
  
  out.preds<-predictSurface.Krig(object = obj4.raw.tt, grid.list = grd,
                                 ZGrid = predgrd,
                                 extrap=TRUE)
  dataOut = c(dataOut, c(out.preds$z))
  x = c(x, rep(out.preds$x, out.preds$ny))
  y = c(y, sort(rep(out.preds$y, out.preds$nx)))
  t = c(t, rep(as.numeric(tsteps[ii]), out.preds$nx*out.preds$ny))
}


AllData = data.frame(preds = dataOut, Lat = y, Lon = x, t = t)
AllData = AllData[!is.na(AllData$preds),]
AllData$LocIDX = paste0(AllData$Lat, AllData$Lon)
# Use the standard Errors to restrict the grid to locations with SE less 
# than 2 dB


SE_Grid = data.frame(se = c(out.SEs$z),
                     Lon =  rep(out.SEs$x, out.SEs$ny),
                     Lat = sort(rep(out.SEs$y, out.SEs$nx)),
                     LocIDX = paste0(sort(rep(out.SEs$y, out.SEs$nx)),
                                     rep(out.SEs$x, out.SEs$ny)))

GridCorrds = data.frame(
  Lon =  rep(out.SEs$x, out.SEs$ny),
  Lat = sort(rep(out.SEs$y, out.SEs$nx)))

AllData$MinDist =NaN

# Only model data points that were wintin 20km of any sensor
for(ii in 1:length(tsteps)){
  
  AllData_idx = which(AllData$t== tsteps[ii])
  
  # Figure out where the sensors are for each of the timesteps
  UTC_timeids  = tseries[1]+ (mean(telapsed)+tsteps[ii]*60)
  
  # Assuming df is your dataframe and given_utc is your single UTC value
  closest_points <- lapply(unique(GPSdf$DriftName), function(id) {
    drift_df <- GPSdf[GPSdf$DriftName == id, ]
    drift_df$diff <- abs(drift_df$UTC - UTC_timeids)
    closest_index <- which.min(drift_df$diff)
    return(drift_df[closest_index, c("DriftName", "Latitude", "Longitude", "UTC")])
  })
  
  # Convert the list to a dataframe
  closest_points_df <- do.call(rbind, closest_points)
  
  # Calculate the distance matrix (in km)
  dist_matrix <- rdist.earth(GridCorrds, 
                             closest_points_df[,c("Longitude", "Latitude")],
                             miles = FALSE)
  
  # Find the minimum distance for each point in AllData
  min_distances <- apply(dist_matrix, 1, min)
  
  # Set NL to NA for points more than 20km away from any drifter
  AllData$MinDist[AllData_idx] <- min_distances
}


# Dummy for preds while working on it
AllData$predsTemp = AllData$preds
AllData$predsTemp[AllData$MinDist>5] = NA
keepIdx = SE_Grid$LocIDX[SE_Grid$se<4]


AllData$Keep = AllData$LocIDX %in% keepIdx

dataSub = AllData[AllData$Keep == TRUE,]
dataSub = AllData[!is.na(AllData$predsTemp),]

dataSub= AllData[AllData$MinDist<5,]

p<-ggplot(subset(dataSub),
          aes(x = Lon, y=Lat, fill = predsTemp))+
  geom_tile()+
  scico::scale_fill_scico(palette = "bilbao", na.value="white")+
  labs(y='Lat', x = 'Lon')+
  theme_bw()


p + transition_time(t) +
  labs(title = "RelativeMinutes: {frame_time}")


anim_save(filename = 'ModelledNL Fields Matern 20km.gif')

############################################################################
# Same thing with the thin plate spline
############################################################################



dataOut = array()
x = array()
y= array()
t = array()

for(ii in 1:length(tsteps)){
  predgrd <- list(x = seq(from = min(loc[,1]), to = max(loc[,1]), length.out = 200),
                  y= seq(from = min(loc[,2]), to = max(loc[,2]), length.out = 400),
                  z = matrix(nrow = 200, ncol = 400, data = tsteps[ii]))
  
  out.preds<-predictSurface(object = obj.TPS.raw, grid.list = grd,
                            ZGrid = predgrd,
                            extrap=TRUE)
  dataOut = c(dataOut, c(out.preds$z))
  x = c(x, rep(out.preds$x, out.preds$ny))
  y = c(y, sort(rep(out.preds$y, out.preds$nx)))
  t = c(t, rep(as.numeric(tsteps[ii]), out.preds$nx*out.preds$ny))
}

AllData = data.frame(preds = dataOut, Lat = y, Lon = x, t = t)
AllData = AllData[!is.na(AllData$preds),]
AllData$LocIDX = paste0(AllData$Lat, AllData$Lon)
# Use the standard Errors to restrict the grid to locations with SE less 
# than 2 dB


SE_Grid = data.frame(se = c(out.SEs$z),
                     Lon =  rep(out.SEs$x, out.SEs$ny),
                     Lat = sort(rep(out.SEs$y, out.SEs$nx)),
                     LocIDX = paste0(sort(rep(out.SEs$y, out.SEs$nx)),
                                     rep(out.SEs$x, out.SEs$ny)))

GridCorrds = data.frame(
  Lon =  rep(out.SEs$x, out.SEs$ny),
  Lat = sort(rep(out.SEs$y, out.SEs$nx)))

AllData$MinDist =NaN

# Only model data points that were wintin 20km of any sensor
for(ii in 1:length(tsteps)){
  
  AllData_idx = which(AllData$t== tsteps[ii])
  
  # Figure out where the sensors are for each of the timesteps
  UTC_timeids  = tseries[1]+ (mean(telapsed)+tsteps[ii]*60)
  
  # Assuming df is your dataframe and given_utc is your single UTC value
  closest_points <- lapply(unique(GPSdf$DriftName), function(id) {
    drift_df <- GPSdf[GPSdf$DriftName == id, ]
    drift_df$diff <- abs(drift_df$UTC - UTC_timeids)
    closest_index <- which.min(drift_df$diff)
    return(drift_df[closest_index, c("DriftName", "Latitude", "Longitude", "UTC")])
  })
  
  # Convert the list to a dataframe
  closest_points_df <- do.call(rbind, closest_points)
  
  # Calculate the distance matrix (in km)
  dist_matrix <- rdist.earth(GridCorrds, 
                             closest_points_df[,c("Longitude", "Latitude")],
                             miles = FALSE)
  
  # Find the minimum distance for each point in AllData
  min_distances <- apply(dist_matrix, 1, min)
  
  # Set NL to NA for points more than 20km away from any drifter
  AllData$MinDist[AllData_idx] <- min_distances
}


# Dummy for preds while working on it
AllData$predsTemp = AllData$preds
AllData$predsTemp[AllData$MinDist>5] = NA

keepIdx = SE_Grid$LocIDX[SE_Grid$se<4]


AllData$Keep = AllData$LocIDX %in% keepIdx

dataSub = AllData[AllData$Keep == TRUE,]
dataSub = AllData[!is.na(AllData$predsTemp),]

dataSub= AllData[AllData$MinDist<5,]

p<-ggplot(subset(dataSub),
          aes(x = Lon, y=Lat, fill = predsTemp))+
  geom_tile()+
  scico::scale_fill_scico(palette = "bilbao", na.value="white")+
  labs(y='Lat', x = 'Lon')+
  theme_bw()



p + transition_time(t) +
  labs(title = "RelativeMinutes: {frame_time}")


p + transition_time(t) +
  labs(title = "RelativeMinutes: {frame_time}")
anim_save(filename = 'ModelledNL Fields ThinPlateSpline.gif')

##########################################################################
# Try GAMs
#########################################################################

library(mgcv)

tseries = simplifNL$UTC[simplifNL$Band== '500Hz'] 
telapsed = (tseries-min(tseries))
# exponential vs Wendland covariance function
tt = as.numeric(telapsed-mean(telapsed))/60

allNl$ElapsedMin = as.numeric(allNl$UTC-mean(allNl$UTC))/60

gam.mod1 = gam(data = allNl[allNl$Band== '500Hz',], NL~ te(Lat, Lon, k=4)+
                s(ElapsedMin, k=3))


gam.mod2 = gam(data = allNl[allNl$Band== '500Hz',], NL~ te(Lat, Lon, k=4)+
                 s(ElapsedMin, k = 5))

gam.mod3 = gam(data = allNl[allNl$Band== '500Hz',], NL~ te(Lat, Lon, k=5)+
                 s(ElapsedMin, k = 5))

gam.mod4 = gam(data = allNl[allNl$Band== '500Hz',], NL~ te(Lat, Lon, k=7)+
                 s(ElapsedMin, k = 7))

gam.mod5 = gam(data = allNl[allNl$Band== '500Hz',], NL~ te(Lat, Lon, k=9)+
                 s(ElapsedMin, k = 9))


# gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ s(Lat, Lon, k= 10)+
#                  s(ElapsedMin, k = 10), family = Gamma(link = "inverse"))
# 
# gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ s(Lat, Lon, k= 10)+
#                  s(ElapsedMin, k = 10), family = Gamma(link = "inverse"))
# 
# 
# 
# 
# gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ te(Lat, Lon, k=10)+
#                  s(ElapsedMin, k = 10), family = Gamma(link = "inverse"))
# 
# gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ te(Lat, Lon, k=5)+
#                  s(ElapsedMin, k = 5), family = Gamma(link = "inverse"))
# 
# gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ te(Lat, Lon, k=7)+
#                  s(ElapsedMin, k = 7), 
#                family = Gamma(link = "inverse"))
# 
# gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ s(Lat, Lon, k=5)+
#                  s(ElapsedMin, k = 5), family = Gamma(link = "inverse"))


gam.mod6 = gam(data = allNl[allNl$Band== '500Hz',], 
               NL~ te(Lat, Lon, k=10)+
                 s(ElapsedMin, k = 350))

# 
# gam.mod7 = gam(data = allNl[allNl$Band== '500Hz',], 
#                NL~ te(Lat, Lon, k=10)+
#                  s(ElapsedMin, k = 300),
#                family = Gamma(link = "inverse"))

AIC(gam.mod1, gam.mod2, gam.mod3, gam.mod4, gam.mod5, gam.mod6,
    gam.mod7)

allNl$GamPred= predict.gam(gam.mod6, newdata = allNl)


#########################################################
# GAM models with fluctutations
#########################################################

# the above method is all well and good but fitting a temporal spline with
# 300 knots seems inefficient at best and stupid overall. One of the tings we
# are actually interested in is variation in space, regardless of the short 
# time period. What if we model the GAM as the difference in mean, I know we
# did that above but wait. 

# Wait no, the model the hourly median noise level as a function of space. Doing
# that we will effectively remove the brief perturbations but it won't sho where
# those perterbations are....

simplifNL$ElapsedMin <- as.numeric((simplifNL$UTC- mean(simplifNL$UTC)))

gam.hrly.01 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
               NL~ s(Lat)+ s(Lon)+
                 s(ElapsedMin, k = 20),
               family = Gamma(link = "inverse"))

gam.hrly.02 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=3)+ s(Lon, k=3)+
                    s(ElapsedMin, k = 20),
                  family = Gamma(link = "inverse"))


gam.hrly.03 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=13)+ s(Lon, k=13)+
                    s(ElapsedMin, k = 20),
                  family = Gamma(link = "inverse"))


gam.hrly.04 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=13)+ s(Lon, k=13)+
                    s(ElapsedMin, k = 10),
                  family = Gamma(link = "inverse"))

gam.hrly.05 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=13)+ s(Lon, k=13)+
                    s(ElapsedMin, k = 30),
                  family = Gamma(link = "inverse"))


gam.hrly.05 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=7)+ s(Lon, k=7)+
                    s(ElapsedMin, k = 40),
                  family = Gamma(link = "inverse"))


gam.hrly.06 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=7)+ s(Lon, k=7)+
                    s(ElapsedMin, k = 50),
                  family = Gamma(link = "inverse"))

gam.hrly.07 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat, k=9)+ s(Lon, k=7)+
                    s(ElapsedMin, k = 50),
                  family = Gamma(link = "inverse"))


gam.hrly.08 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat,Lon, k=7)+
                    s(ElapsedMin, k = 50),
                  family = Gamma(link = "inverse"))


gam.hrly.09 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat,Lon, k=9)+
                    s(ElapsedMin, k = 50),
                  family = Gamma(link = "inverse"))


gam.hrly.10 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat,Lon, k=10)+
                    s(ElapsedMin, k = 50),
                  family = Gamma(link = "inverse"))


gam.hrly.11 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat,Lon, k=15)+
                    s(ElapsedMin, k = 50),
                  family = Gamma(link = "inverse"))


gam.hrly.12 = gam(data = simplifNL[simplifNL$Band== '500Hz',], 
                  NL~ s(Lat,Lon, ElapsedMin, k=120),
                  family = Gamma(link = "inverse"))


gam.hrly.13 <- gam(NL ~ s(Lat, Lon,by = ElapsedMin), 
                 data = simplifNL[simplifNL$Band== '500Hz',])

gam.hrly.14 <- gam(NL ~ te(Lat, Lon) + s(ElapsedMin,k=120),
                   data = simplifNL[simplifNL$Band== '500Hz',],
                   family = Gamma(link = "inverse"))

gam.hrly.14 <- gam(NL ~ te(Lat, Lon, ElapsedMin),
                   data = simplifNL[simplifNL$Band== '500Hz',])

gam.hrly.15 <- gam(RelNoise ~ te(Lat, Lon),
                   data = simplifNL[simplifNL$Band== '500Hz',])


AIC(gam.hrly.01, gam.hrly.02, gam.hrly.03, gam.hrly.04, gam.hrly.05, 
    gam.hrly.06, gam.hrly.07, gam.hrly.08, gam.hrly.09, gam.hrly.10,
    gam.hrly.11, gam.hrly.12,gam.hrly.13, gam.hrly.14, gam.hrly.15)


plot(gam.hrly.15, resid = TRUE)



# Create the prediction dataframe and clean it out 
PredDF = expand.grid(
  Lon = out.SEs$x, 
  Lat = out.SEs$y)

# Determine the minimum approc distance for each grid in order to make a mask
dist_matrix <- rdist.earth(as.matrix(PredDF[, c("Lon", "Lat")]), 
                           as.matrix(GPSdf[, c("Longitude", "Latitude")]), miles = FALSE)
# Grid locations more than 10km from a gps ping
PredDF$Mask <- apply(dist_matrix, 1, function(x) any(x <= 10))

PredDF = PredDF[PredDF$Mask==TRUE,]

timeStamps = sort(
  rep(
    quantile(simplifNL$ElapsedMin, seq(0,1, by=.02)),
    nrow(PredDF)))

PredDF= PredDF[rep(1:nrow(PredDF),51),]
PredDF$ElapsedMin = timeStamps

# Make the predictions
PredDF$NLpred = predict(gam.hrly.15, PredDF,type = "response" )

# Add the GPS elapsed time
GPSdf$ElapsedMin <- as.numeric((GPSdf$UTC- mean(simplifNL$UTC)))


# Create a new GPS dataframe that is binned by the timestamps
tsteps = unique(PredDF$ElapsedMin)




datasub = subset(GPSdf, DriftName == DriftNames[1])
df = data.frame(
  ElapsedMin= tsteps,
  Lat = approx(datasub$ElapsedMin, datasub$Latitude, tsteps)$y,
  Lon = approx(datasub$ElapsedMin, datasub$Longitude, tsteps)$y)



for(ii in 2:length(DriftNames)){
  
  datasub = subset(GPSdf, DriftName == DriftNames[ii])
  df = rbind(df, data.frame(
    ElapsedMin= tsteps,
    Lat = approx(datasub$ElapsedMin, datasub$Latitude, tsteps)$y,
    Lon = approx(datasub$ElapsedMin, datasub$Longitude, tsteps)$y))
  
}

# Make an initial plot
p<-ggplot()+
  geom_tile(data = PredDF,
            aes(x = Lon, y=Lat, fill = NLpred))+
  scico::scale_fill_scico(palette = "bilbao", na.value="white")+
  geom_point(data = GPSdf[,1:16], aes(x=Longitude, y=Latitude), color='black')+
  #geom_point(data = df, aes(x=Lon, y=Lat), color='red', size =2)+
  labs(y='Latitude', x = 'Lonitude')+
  theme_bw()


p + transition_time(ElapsedMin) +
  labs(title = "RelativeMinutes: {frame_time}")

anim_save(filename = 'ModelledNL Fields latxlonxtime.gif')

 

###############################################
# Plot the predictions
###############################################


tsteps = quantile(allNl$ElapsedMin, seq(0,1, by=.2))

allData = data.frame(
  Lon =  rep(rep(out.SEs$x, out.SEs$ny), length(tsteps)),
  Lat = rep(sort(rep(out.SEs$y, out.SEs$nx))), length(tsteps),
  ElapsedMin = sort(rep(tsteps, out.SEs$nx*out.SEs$ny)))

allData$MinDist= NA
allData$NLpred = NA



# Only model data points that were wintin 20km of any sensor
for(ii in 1:length(tsteps)){
  
  AllData_idx = which(allData$ElapsedMin== tsteps[ii])
  
  # Figure out where the sensors are for each of the timesteps
  UTC_timeids  = mean(allNl$UTC)+ (tsteps[ii]*60)
  
  # Assuming df is your dataframe and given_utc is your single UTC value
  closest_points <- lapply(unique(GPSdf$DriftName), function(id) {
    drift_df <- GPSdf[GPSdf$DriftName == id, ]
    drift_df$diff <- abs(drift_df$UTC - UTC_timeids)
    closest_index <- which.min(drift_df$diff)
    return(drift_df[closest_index, c("DriftName", "Latitude", "Longitude", "UTC")])
  })
  
  # Convert the list to a dataframe
  closest_points_df <- do.call(rbind, closest_points)
  
  # Calculate the distance matrix (in km)
  dist_matrix <- rdist.earth(GridCorrds, 
                             closest_points_df[,c("Longitude", "Latitude")],
                             miles = FALSE)
  
  # Find the minimum distance for each point in AllData
  min_distances <- apply(dist_matrix, 1, min)
  
  # Set NL to NA for points more than 20km away from any drifter
  allData$MinDist[AllData_idx] <- min_distances
  
  allData$NLpred[AllData_idx]= predict.gam(gam.mod6, 
                                         newdata =  
                                           allData[AllData_idx,],
                                         type='response')

}

ggplot(subset(allData, MinDist<=15),
       aes(y= Lat, x= Lon))+
  facet_wrap(~ElapsedMin, nrow=3)+
  geom_tile(aes(fill= NLpred))

ggplot(subset(allData, MinDist<=10 & NLpred>65 & NLpred<130),
       aes(y= Lat, x= Lon))+
  facet_wrap(~ElapsedMin, nrow=3)+
  geom_tile(aes(fill= NLpred))


# Dummy for preds while working on it
AllData$predsTemp = AllData$preds
AllData$predsTemp[AllData$MinDist>5] = NA

keepIdx = SE_Grid$LocIDX[SE_Grid$se<4]


AllData$Keep = AllData$LocIDX %in% keepIdx

dataSub = AllData[AllData$Keep == TRUE,]
dataSub = AllData[!is.na(AllData$predsTemp),]

dataSub= AllData[AllData$MinDist<5,]

p<-ggplot(subset(dataSub),
          aes(x = Lon, y=Lat, fill = predsTemp))+
  geom_tile()+
  scico::scale_fill_scico(palette = "bilbao", na.value="white")+
  labs(y='Lat', x = 'Lon')+
  theme_bw()


#########################################################################
# Gstats has spatial-temporal fitting
##########################################################################
library(gstat)
library(spacetime)
library(sp)

# Try to convert the data to a spacetime dataframe

spNoise = SpatialPoints(simplifNL[simplifNL$Band=='500Hz',c()])

noiseSTDF



rn = row.names(noise500reshape)[2:7]
acf(noise500reshape[,2:8],  na.action = na.pass)
acf(noise20kreshape[,2:8],  na.action = na.pass)

acf(na.omit(as(r5to10[rn,], "xts")))


vv = variogram(PM10~1, r5to10, width=20, cutoff = 200, tlags=0:5)






































