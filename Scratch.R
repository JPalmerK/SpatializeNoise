rm(list = ls())
library(mgcv)
library(dplyr)
library(lubridate)
library(viridis)
library(ggplot2)
library(tidyr)
library(ggcorrplot)
library(sp)

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




#load('AllNL')


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
obj.lf.norm<- spatialProcess( loc, nl.norm.lf)
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
out.exponential <- optim(par=c(2,50,1),fn=ls.exponential,method="L-BFGS", lower = c(0,0,0), upper = c(80,100,80))

rho.exponential <- out.exponential$par[1]
aRange.exponential <- out.exponential$par[2]
sigma2.exponential <- out.exponential$par[3]

plot(v.lf$stats["mean",1:17]~v.lf$centers[1:17], main="Exponential vs. Matern")
lines(c(sigma2.exponential + rho.exponential*(1-exp(-v.lf$centers/aRange.exponential)))~v.lf$centers,col=1)

# Mattern plot isn't working so we will skip it.
###############################################################################
# Using spatila process to fit spatial models, rather than step by step covariance
# selection
######################################################################

fit.lf.norm<- spatialProcess(x = loc, y= nl.norm.lf)
surface(fit.lf.norm)
title("Kriging Predictions")

# plot the fitted model to evaluate fit
set.panel(2,2)
plot(fit.lf.norm)

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


# exponential vs Wendland covariance function
tt = as.numeric(telapsed-mean(telapsed))/60

obj.norm <- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth")


obj2<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
                       cov.args = list(Covariance ="Exponential")) #nope

obj3<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
                       cov.args = list(Covariance = "Wendland",
                                       dimension = 1, k = 1))

obj4<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
                       cov.args = list(Covariance = "Matern"))

rbind(obj.norm$summary, obj2$summary, obj3$summary, obj3$summary)
# Plot 1 shows data vs. predicted values, and plot 2 shows predicted values vs.
# residuals. Plot 3 shows the criteria to select the smoothing parameter λ = σ
# 2/ρ
set.panel(1,2)
plot(obj.norm)
plot(obj2)
plot(obj3)
plot(obj4)


obj.raw.tt <- spatialProcess( x = loc, y = nl.lf, Z = tt,
                              Distance = "rdist.earth")
obj2.raw.tt<- spatialProcess( loc, nl.lf, Distance = "rdist.earth",
                       cov.args = list(Covariance ="Exponential"))#nope

obj3.raw.tt<- spatialProcess( loc,  y = nl.lf, Z = tt,
                              Distance = "rdist.earth",
                       cov.args = list(Covariance = "Wendland",
                                       dimension = 1, k = 1))

obj4.raw.tt<- spatialProcess( loc, y = nl.lf, Z = tt,
                              Distance = "rdist.earth",
                       cov.args = list(Covariance = "Matern"))


rbind(obj.raw.tt$summary,obj2.raw.tt$summary,obj3.raw.tt$summary,
      obj4.raw.tt$summary)

set.panel(1,2)
plot(obj.raw.tt)
plot(obj2.raw.tt)
plot(obj3.raw.tt)
plot(obj4.raw.tt)#best


# Samething but normalized
obj.norm.tt <- spatialProcess( x = loc, y = nl.norm.lf, Z = tt,
                              Distance = "rdist.earth")
obj2.norm.tt<- spatialProcess( loc, nl.norm.lf, Distance = "rdist.earth",
                              cov.args = list(Covariance ="Exponential"))

obj3.norm.tt<- spatialProcess( loc,  y = nl.norm.lf, Z = tt,
                              Distance = "rdist.earth",
                              cov.args = list(Covariance = "Wendland",
                                              dimension = 1, k = 1))

obj4.norm.tt<- spatialProcess( loc, y = nl.norm.lf, Z = tt,
                              Distance = "rdist.earth",
                              cov.args = list(Covariance = "Matern"))


rbind(obj.norm.tt$summary,
      obj2.norm.tt$summary,
      obj3.norm.tt$summary,
      obj4.norm.tt$summary)

set.panel(1,2)
plot(obj.norm.tt)
plot(obj2.norm.tt)
plot(obj3.norm.tt)
plot(obj4.norm.tt)


###########################################################################
# Plot the best model
###########################################################################



obj4.raw.tt<- spatialProcess( loc, y = nl.lf, Z = tt,
                              Distance = "rdist.earth",
                              cov.args = list(Covariance = "Matern"))

predgrd <- list(x = seq(from = min(loc[,1]), to = max(loc[,1]), length.out = 200),
            y= seq(from = min(loc[,2]), to = max(loc[,2]), length.out = 400),
            z = matrix(nrow = 200, ncol = 400, data = -3660.63981))


out.preds<-predictSurface.Krig(object = obj4.raw.tt, grid.list = grd,
                            ZGrid = predgrd, 
                            extrap=TRUE) 

image.plot( out.preds, col=larry.colors())

out.SEs<-predictSurfaceSE(object = obj4.raw.tt, grid.list = grd,
                            ZGrid = predgrd, 
                            extrap=TRUE) 


