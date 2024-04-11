# MRSEA modelling
rm(list = ls())
library(tidyr)
library(dplyr)
library(MRSea)
library(sp)
library(ggplot2)
library(lubridate)

#\source('//getRadiiSequence.R')
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



allNl <- allNl %>%
  mutate(Date = as.Date(UTC),
         Hour = hour(UTC))

allNl = allNl[!is.na(allNl$Lat),]

# Trim data where instruments were pulled
allNl=allNl[allNl$Lat>35.2,]


# Get UTM values
library(sp)
cord.dec = SpatialPoints(cbind(allNl$Lon, allNl$Lat),
                         proj4string=CRS("+proj=longlat"))

out = spTransform(cord.dec, CRS("+init=epsg:32610"))

allNl$x.pos = coordinates(out)[,1]
allNl$y.pos = coordinates(out)[,2]

# Calculate the median of 'NL' for each combination of 'Date' and 'Hour'
allNl <- allNl %>%
  group_by(Band, Date, Hour) %>%
  mutate(Median_NL = median(NL))


# Detrend the data
allNl <- allNl %>%
  group_by(Band) %>%
  dplyr::arrange(desc(UTC)) %>% 
  dplyr::mutate(NL_ave = zoo::rollmedian(NL, k = 211, fill = NA)) %>%
  dplyr::mutate(NL_ave_mean = 20*log10(zoo::rollmean(10^(NL/20), k = 211, fill = NA))) %>%
  ungroup()


# Model response
allNl$response = allNl$NL- allNl$NL_ave_mean
allNl$Band <- factor(allNl$Band , levels = c("20kHz", "500Hz"))

# Plot data and trends       
ggplot(allNl)+
  geom_line(aes(x=UTC, y=NL, color = DriftName))+
  scale_color_brewer(palette = 'Paired', name= 'Drift ID')+
  geom_line(aes(x=UTC, y= NL_ave_mean), color ='black')+
  facet_wrap(~Band, nrow = 2)+
  theme_bw()+ylab('Third-Octave NL')+xlab('UTC')


lfNL = subset(allNl, Band == '500Hz')
hfNL = subset(allNl, Band == '20kHz')


ggplot(hfNL[idx_critters,])+
  geom_line(aes(x=UTC, y=NL, color = DriftName))+
  scale_color_brewer(palette = 'Paired', name= 'Drift ID')+
  geom_line(aes(x=UTC, y= NL_ave_mean), color ='black')+
  theme_bw()+ylab('Third-Octave NL')+xlab('UTC')

rm(list=setdiff(ls(), c("allNl", 'lfNL', 'hfNL')))

###################################################################
# Look at standard deviation between instruments and between the baseline
##################################################################

# If the spatial component is important, we should see a considerable 
# variation between units. This should be greater than the variation in time

# Detrend the data
allNl <- allNl %>%
  group_by(Band) %>%
  dplyr::arrange(desc(UTC)) %>% 
  dplyr::mutate(NL_ave = zoo::rollv(NL, k = 211, fill = NA)) %>%
  dplyr::mutate(NL_ave_mean = 20*log10(zoo::rollmean(10^(NL/20), k = 211, fill = NA))) %>%
  ungroup()






###################################################################
# Check for co-linearity
##################################################################

hfNL$seconds = as.numeric(hfNL$UTC-median(hfNL$UTC))

covariates <- c("x.pos", "y.pos", "seconds","response")
pairs(subset(hfNL, select=covariates), 
      upper.panel=NULL, pch=19, cex=0.3)


hfNL = hfNL[!is.na(hfNL$response),]

hfGLM <- glm(response ~ x.pos + y.pos+seconds ,
             data=hfNL)
car::vif(hfGLM)

###############################################################
# Fit an initial model and look at residuals
###############################################################
hfNL = hfNL[!is.na(hfNL$response),]

hfGLM <- glm(response ~ 1 , data=hfNL)

mdlSummary <- data.frame(Observed=hfGLM$model$response,
                         Fitted=predict(hfGLM,  type="response"), 
                         Residuals=residuals(hfGLM, type="pearson"),
                         Index=seq(length(hfGLM$model$response)))
ggplot(mdlSummary) +
  geom_line(aes(x=Index, y=Observed, col="Observed"), lwd=1) +
  geom_line(aes(x=Index, y=Fitted, col="Fitted"), lwd=1) +
  scale_color_manual(values=c('Observed'="#377eb8", 'Fitted'="#4daf4a")) +
  labs(color="") 

ggplot(mdlSummary) +
  geom_line(aes(x=Index, y=Residuals, col="Residuals"), lwd=1) +
  scale_color_manual(values=c('Residuals'="#e41a1c")) +
  labs(color="")  +
  ylab("Bird counts")


# Perform Fourier analysis to see if there is a tidal cycle
fftdata =hfNL[hfNL$DriftName=='ADRIFT_046', ]
fft_result <- fft(fftdata$NL_ave)
freq <- seq(0, 1/(2 * 60*as.numeric(fftdata$UTC[2] - fftdata$UTC[1])), 
            length.out = length(fft_result))


# Plot the amplitude spectrum
plot(freq[1:800], log10(Mod(fft_result[1:800])),
     type = "l", xlab = "Frequency", 
     ylab = "Amplitude", main = "Fourier Transform")


##########################################################################
# Set the knot grid
##########################################################################
set.seed(123)

knotgrid<- getKnotgrid(coordData = cbind(hfNL$Lon, hfNL$Lat),
                         numKnots = 450,
                         plot = TRUE)

distMatshavers <- makeDistsHaversine(cbind(hfNL$Lon, hfNL$Lat), knotgrid)
dimnames(distMatshavers$knotDist)[[1]]<-as.character(1:450)
dimnames(distMatshavers$knotDist)[[2]]<-as.character(1:450)


distMats <- makeDistsHaversine(cbind(hfNL$Lon, hfNL$Lat), knotgrid)

# Convert the corrdinates to UTM
cord.dec = SpatialPoints(cbind(knotgrid[,1], knotgrid[,2]), 
                         proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32610"))
knotgridhavers_converted =as.matrix(coordinates(cord.UTM))
colnames(knotgridhavers_converted)





################################################################
#Salsa model
################################################################

salsa2dlist<-list(fitnessMeasure = 'BIC',
                  knotgrid = knotgridhavers_converted,
                  minKnots_1d = c(5), 
                  maxKnots_1d = c(40), 
                  startKnots_1d = c(20), 
                  degree = c(NA),
                  gaps = c(0),
                  removal = TRUE,
                  splines = c("ns"))
# Gaussian
salsa2dOutput<-runSALSA2D(model = hfGLM,
                          salsa2dlist = salsa2dlist,
                          d2k=distMats$dataDist,
                          k2k=distMats$knotDist,
                          suppress.printout = FALSE)
# Exponential basis
salsa2dOutput.exp<-runSALSA2D(hfGLM,
                              salsa2dlist, 
                              d2k=distMats$dataDist,
                              k2k=distMats$knotDist,
                              basis = "exponential", ##
                              suppress.printout = FALSE)


# Gaussian log transformed distance because I'm a chaos demon!
salsa2dOutput_log<-runSALSA2D(model = hfGLM,
                          salsa2dlist = salsa2dlist,
                          d2k=log10(distMats$dataDist),
                          k2k=log10(distMats$knotDist),
                          suppress.printout = FALSE)



ggplot(hfNL) + 
  geom_tile(aes(x=x.pos, y=y.pos, 
                fill=fitted(salsa2dOutput$bestModel), 
                height=300, width=300)) +
  scale_fill_distiller(palette = "Spectral",name="Noise Level") +
  xlab("Easting (km)") + ylab("Northing (km)") + 
  theme_bw()


ggplot(hfNL) + 
  geom_tile(aes(x=x.pos, y=y.pos, 
                fill=fitted(salsa2dOutput.exp$bestModel), 
                height=300, width=300)) +
  scale_fill_distiller(palette = "Spectral",name="Noise Level") +
  xlab("Easting (km)") + ylab("Northing (km)") + 
  theme_bw()


ggplot(hfNL) + 
  geom_tile(aes(x=x.pos, y=y.pos, 
                fill=fitted(salsa2dOutput_log$bestModel), 
                height=300, width=300)) +
  scale_fill_distiller(palette = "Spectral",name="Noise Level") +
  xlab("Easting (km)") + ylab("Northing (km)") + 
  theme_bw()



# Raw data
ggplot(hfNL[15043:17500,]) + 
  geom_point(aes(x=Lon, y=Lat, color=NL)) +
  scale_color_distiller(palette = "Spectral",name="Noise Level") +
  xlab("Easting (km)") + ylab("Northing (km)") + 
  theme_bw()

####################################################
# Make Prediction Grid
####################################################
# Prediction Space
nlPredData = expand.grid(Lon = 
                                  seq(min(allNl$Lon), max(allNl$Lon), 
                                      length.out =200),
                                Lat = 
                                  seq(min(allNl$Lat), max(allNl$Lat), 
                                      length.out =500))
# Convert the corrdinates to UTM
cord.dec = SpatialPoints(cbind(nlPredData$Lon, nlPredData$Lat), 
                         proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32610"))

nlPredData$x.pos = coordinates(cord.UTM)[,1]
nlPredData$y.pos = coordinates(cord.UTM)[,2]


# Create the distance mat for the prediction grid
preddist <-makeDistsHaversine(
  cbind(nlPredData$Lon, nlPredData$Lat),
  knotgrid, knotmat=FALSE)$dataDist



# make predictions on response scale
preds<-predict(newdata = nlPredData,
               g2k = preddist,
               object = salsa2dOutput$bestModel)
nlPredData$preds<-preds[,1]
nlPredData$response = nlPredData$preds


nlPredData$predsexp<-predict(newdata = nlPredData,
                             g2k = preddist_havers,
                             object = salsa2dOutput.exp$bestModel)[,1]
nlPredData$mindist =(apply(X = preddist, 1, min))

#Chaaaaaaaaaaaos deamon
nlPredData_log=nlPredData # No this isn't right
nlPredData$preds_log<-predict(newdata = nlPredData,
                             g2k = log10(preddist_havers),
                             object = salsa2dOutput_log$bestModel)[,1]



idx_critters = which(hfNL$UTC> as.POSIXct('2023-03-15 00:00:00', tz= 'UTC') &
                       hfNL$UTC< as.POSIXct('2023-03-15 06:00:00', tz= 'UTC'))
predIDX = range(hfNL$Lat[idx_critters])


idx_prestorm = which(hfNL$UTC< as.POSIXct('2023-03-14 00:00:00', tz= 'UTC'))
predIDX = range(hfNL$Lat[idx_prestorm])

ggplot(nlPredData[nlPredData$mindist<2300 & 
                    nlPredData$Lat>predIDX[1] &  
                   nlPredData$Lat<predIDX[2],]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=preds)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle('SALSA Model Predictions Gaussian Basis')+
  theme(plot.title = element_text(hjust = 0.5))+ theme_bw()+
  geom_point(data = hfNL[idx_critters,],aes(x=Lon, y=Lat, color= NL)) +
  #scale_color_distiller(palette = "Spectral",name="Noise Level") +
  xlab("Easting (km)") + ylab("Northing (km)") 
 
ggplot(hfNL)+
  geom_point(data = hfNL[idx_critters,],aes(x=Lon, y=Lat, color= response)) +
  scale_color_distiller(palette = "Spectral",name="Noise Level") +
  xlab("Easting (km)") + ylab("Northing (km)") 

ggplot(nlPredData[nlPredData$mindist<5000,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=predsexp)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()



ggplot(nlPredData[nlPredData$mindist<5000,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=preds_log)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()





################################################################
# Make a list
################################################################

endTime = min(allNl$UTC[allNl$Lat>35])
beginTime = as.POSIXct('2023-03-11 22:20:00')
endTime = as.POSIXct('"2023-03-11 19:00:00 UTC"')


################################################################
# Use variogram to get the central tendency
################################################################

#https://github.com/lindesaysh/MRSea/blob/HEAD/R/getRadiiSequence.R



rs<-getRadiiSequence(method = "variogram",
                     numberofradii = 8,
                     xydata = hfNL[,c("x.pos", "y.pos")],
                     response = hfNL$response,
                     basis = "gaussian",
                     distMatrix = distMats$dataDist)

rs.exp<-getRadiiSequence(method = "variogram",
                         numberofradii = 8,
                         xydata = hfNL[,c("x.pos", "y.pos")],
                         response = hfNL$response,
                         basis = "exponential",
                         distMatrix = distMats$dataDist)




# In this case, the range parameter is 10.09 which suggests on average, 
# the spatial correlation decays after approximately 4.5km. Which is equivalent
# to 18*log10(4567) or 65 dB


salsa2dlist$r_seq <- rs

salsa2dOutput.vario <-runSALSA2D(model = hfGLM,
                                 salsa2dlist = salsa2dlist,
                                 d2k=distMats$dataDist,
                                 k2k=distMats$knotDist,
                                 basis = 'gaussian',
                                 suppress.printout = FALSE)



# Exponentaial
salsa2dlist$r_seq <- rs.exp
salsa2dOutput.vario.exp <-runSALSA2D(model = hfGLM,
                                     salsa2dlist = salsa2dlist,
                                     d2k=distMats$dataDist,
                                     k2k=distMats$knotDist,
                                     basis = "exponential", ##
                                     suppress.printout = FALSE)


####################################################################
# Assess Model Fits
##################################################################

cv.gamMRSea(hfNL, salsa2dOutput.vario$bestModel, K=10)$delta[2]
cv.gamMRSea(hfNL, salsa2dOutput$bestModel, K=25)$delta[2]

cv.gamMRSea(hfNL, salsa2dOutput.exp$bestModel, K=10)$delta[2]
cv.gamMRSea(hfNL, salsa2dOutput.vario.exp$bestModel, K=10)$delta[2]


# The CV is lower on the exponential but it doesn't look right and 
# it chose more knots. Even thought the CV was higher


hfGLM_nl <- glm(response ~ 1+offset(NL_ave), data=hfNL)


# Plot Pearson residuals vs covariate
plot(hfNL$x.pos, residuals(hfGLM_nl, type="pearson"),
     ylab="Pearson residuals", xlab="Covariate (x)", pch=19)
abline(h=0, lty=2, col="lightgrey", lwd=4)

###############################################################
# Plot the variogram fits
##############################################################


nlPredData$predsexpvario<-predict(newdata = nlPredData,
                                  g2k = preddist,
                                  object = salsa2dOutput.vario.exp$bestModel)


nlPredData$predsvario<-predict(newdata = nlPredData,
                               g2k = preddist,
                               object = salsa2dOutput.vario$bestModel)[,1]



ggplot(nlPredData[nlPredData$mindist<5000,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=predsexpvario)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()


ggplot(nlPredData[nlPredData$mindist<5000,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=predsvario)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()

##########################################################


