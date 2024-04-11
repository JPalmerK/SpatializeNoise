# MRSEA modelling
rm(list = ls())
library(tidyr)
library(dplyr)
library(MRSea)
library(sp)
library(ggplot2)
library(lubridate)

#source('//getRadiiSequence.R')
########################################################################
# Load data, drift GPS, and noise levels pre-process, and wind lease area
########################################################################
# GPS Directory
# csv_directory <- "F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023"
csv_directory <- "C:/Users/kaitlin.palmer/Documents/GitHub/NL_spatialize/SpatializeNoise/GPS_CSV-20230923T045356Z-001/MorroBay Mar 2023/"

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
  ungroup()
           

# Model response
allNl$response = allNl$NL- allNl$NL_ave


# Plot data and trends       
ggplot(allNl)+
  geom_line(aes(x=UTC, y=NL, color = DriftName))+
  scale_color_brewer(palette = 'Paired', name= 'Drift ID')+
  geom_line(aes(x=UTC, y= NL_ave), color ='black')+
  facet_wrap(~Band, nrow = 2)+
  theme_bw()



ggplot(allNl)+
  geom_line(aes(x=UTC, y=NL, color = response))+
  scale_color_brewer(palette = 'Paired', name= 'Drift ID')+
  facet_wrap(~Band, nrow = 2)+
  theme_bw()



lfNL = subset(allNl, Band == '500Hz')
hfNL = subset(allNl, Band == '20kHz')



rm(list=setdiff(ls(), c("allNl", 'lfNL', 'hfNL')))
###############################################################
# Fit an initial model and look at residuals
###############################################################
hfNL = hfNL[!is.na(hfNL$response),]

hfGLM <- glm(response ~ 1 ,
             data=hfNL)

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

lsknotgrid<- getKnotgrid(coordData = cbind(hfNL$x.pos, hfNL$y.pos),
                         numKnots = 450,
                         plot = TRUE)


lfdistMats <- makeDists(cbind(hfNL$x.pos, hfNL$y.pos), lsknotgrid)

################################################################
#Salsa model
################################################################

salsa2dlist<-list(fitnessMeasure = 'BIC',
                  knotgrid = lsknotgrid,
                  startKnots=20,
                  minKnots=10,
                  maxKnots=40,
                  gap=0, 
                  splines = c("ns"))

salsa2dOutput<-runSALSA2D(model = hfGLM,
                          salsa2dlist = salsa2dlist,
                          d2k=lfdistMats$dataDist,
                          k2k=lfdistMats$knotDist,
                          suppress.printout = FALSE)
# Exponential basis
salsa2dOutput.exp<-runSALSA2D(hfGLM,
                              salsa2dlist, 
                              d2k=lfdistMats$dataDist,
                              k2k=lfdistMats$knotDist,
                              basis = "exponential", ##
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




salsa2dlist<-list(fitnessMeasure = 'BIC',
                  knotgrid = lsknotgrid,
                  startKnots=20,
                  minKnots=10,
                  maxKnots=40,
                  gap=0, 
                  splines = c("ns"))

salsa2dOutput<-runSALSA2D(model = hfGLM,
                          salsa2dlist = salsa2dlist,
                          d2k=lfdistMats$dataDist,
                          k2k=lfdistMats$knotDist,
                          suppress.printout = FALSE)


####################################################
# Make Prediction Grid
####################################################

nlPredData = expand.grid(x.pos = 
                           seq(min(allNl$x.pos), max(allNl$x.pos), 
                               length.out =200),
                         y.pos = 
                           seq(min(allNl$y.pos), max(allNl$y.pos), 
                               length.out =500))

preddist<-makeDists(cbind(nlPredData$x.pos, nlPredData$y.pos),
                    lsknotgrid, knotmat=FALSE)$dataDist




# make predictions on response scale
preds<-predict(newdata = nlPredData,
               g2k = preddist,
               object = salsa2dOutput$bestModel)

nlPredData$preds<-preds
nlPredData$predsexp<-predict(newdata = nlPredData,
                             g2k = preddist,
                             object = salsa2dOutput.exp$bestModel)

nlPredData$mindist =(apply(X = preddist, 1, min))
nlPredData$response = nlPredData$preds

ggplot(nlPredData[nlPredData$mindist<2300,]) + 
  geom_tile(aes(x=x.pos, y=y.pos, fill=response)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()+
  geom_point(data = hfNL, aes(x=x.pos, y=y.pos, color = response))+
  scale_color_distiller(palette = "Spectral",name="NL Variation")+
  ggtitle('Gaussian Basis')


ggplot(nlPredData[nlPredData$mindist<2300,]) + 
  geom_tile(aes(x=x.pos, y=y.pos, fill= predsexp)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  # geom_point(data = hfNL, aes(x=x.pos, y=y.pos, color = response))+
  # scale_color_distiller(palette = "Spectral",name="NL Variation")+
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()+
  ggtitle('Exponential Basis')


################################################################
# Use variogram to get the central tendency
################################################################

#https://github.com/lindesaysh/MRSea/blob/HEAD/R/getRadiiSequence.R



rs<-getRadiiSequence(method = "variogram",
                     numberofradii = 8,
                     xydata = hfNL[,c("x.pos", "y.pos")],
                     response = hfNL$response,
                     basis = "gaussian",
                     distMatrix = lfdistMats$dataDist)

rs.exp<-getRadiiSequence(method = "variogram",
                         numberofradii = 8,
                         xydata = hfNL[,c("x.pos", "y.pos")],
                         response = hfNL$response,
                         basis = "exponential",
                         distMatrix = lfdistMats$dataDist)

# In this case, the range parameter is 10.09 which suggests on average, 
# the spatial correlation decays after approximately 4.5km. Which is equivalent
# to 18*log10(4567) or 65 dB


salsa2dlist$r_seq <- rs

salsa2dOutput.vario <-runSALSA2D(model = hfGLM,
                                 salsa2dlist = salsa2dlist,
                                 d2k=lfdistMats$dataDist,
                                 k2k=lfdistMats$knotDist,
                                 basis = 'gaussian',
                                 suppress.printout = FALSE)



# Exponentaial
salsa2dlist$r_seq <- rs.exp
salsa2dOutput.vario.exp <-runSALSA2D(model = hfGLM,
                                     salsa2dlist = salsa2dlist,
                                     d2k=lfdistMats$dataDist,
                                     k2k=lfdistMats$knotDist,
                                     basis = "exponential", ##
                                     suppress.printout = FALSE)


####################################################################
# Assess Model Fits
##################################################################

cv.gamMRSea(hfNL, salsa2dOutput.vario$bestModel, K=10)$delta[2]
cv.gamMRSea(hfNL, salsa2dOutput$bestModel, K=10)$delta[2]

cv.gamMRSea(hfNL, salsa2dOutput.exp$bestModel, K=10)$delta[2]
cv.gamMRSea(hfNL, salsa2dOutput.vario.exp$bestModel, K=10)$delta[2]




###############################################################
# Plot the variogram fits
##############################################################


nlPredData$predsexpvario<-predict(newdata = nlPredData,
                                  g2k = preddist,
                                  object = salsa2dOutput.vario.exp$bestModel)


nlPredData$predsvario<-predict(newdata = nlPredData,
                                  g2k = preddist,
                                  object = salsa2dOutput.vario$bestModel)


ggplot(nlPredData[nlPredData$mindist<5000,]) + 
  geom_tile(aes(x=x.pos, y=y.pos, fill=predsexpvario)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()


ggplot(nlPredData[nlPredData$mindist<2500,]) + 
  geom_tile(aes(x=x.pos, y=y.pos, fill=(predsvario))) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()



################################################################
# Try it with geodesic fitting-  working
################################################################

# Knots and distance grids
knotgridhavers<- getKnotgrid(coordData = cbind(hfNL$Lon, hfNL$Lat),
                         numKnots = 250,
                         plot = TRUE)


distMatshavers <- makeDistsHaversine(cbind(hfNL$Lon, hfNL$Lat), knotgridhavers)
dimnames(distMatshavers$knotDist)[[1]]<-as.character(1:250)
dimnames(distMatshavers$knotDist)[[2]]<-as.character(1:250)


# Convert the corrdinates to UTM
cord.dec = SpatialPoints(cbind(knotgridhavers[,1], knotgridhavers[,2]), 
                         proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32610"))

head(knotgridhavers)

knotgridhavers_converted =as.matrix(coordinates(cord.UTM))
colnames(knotgridhavers_converted)


salsa2dlist<-list(fitnessMeasure = 'BIC',
                  knotgrid = knotgridhavers_converted,
                  startKnots=20,
                  minKnots=10,
                  maxKnots=40,
                  gap=0, 
                  splines = c("ns"))


# Gaussian basis
salsa2dOutput.haversine<-runSALSA2D(hfGLM,
                          salsa2dlist, 
                          d2k=distMatshavers$dataDist,
                          k2k=distMatshavers$knotDist,
                          basis = "gaussian", ##
                          suppress.printout = FALSE)




####################################################
# Make Prediction Grid
####################################################

# Prediction Space
nlPredData_havers = expand.grid(Lon = 
                                  seq(min(allNl$Lon), max(allNl$Lon), 
                                      length.out =200),
                                Lat = 
                                  seq(min(allNl$Lat), max(allNl$Lat), 
                                      length.out =500))
# Convert the corrdinates to UTM
cord.dec = SpatialPoints(cbind(nlPredData_havers$Lon, nlPredData_havers$Lat), 
                         proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32610"))

nlPredData_havers$x.pos = coordinates(cord.UTM)[,1]
nlPredData_havers$y.pos = coordinates(cord.UTM)[,2]


# Create the distance mat for the prediction grid
preddist_havers <-makeDistsHaversine(
  cbind(nlPredData_havers$Lon, nlPredData_havers$Lat),
  knotgridhavers, knotmat=FALSE)$dataDist



# make predictions on response scale
preds<-predict(newdata = nlPredData_havers,
               g2k = preddist_havers,
               object = salsa2dOutput$bestModel)




nlPredData_havers$preds <-preds[,1]
ggplot(nlPredData_havers[nlPredData$mindist<5000,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=preds)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()



#########################################################################
# Exponential response
#########################################################################

# Exponential basis
salsa2dOutput.havers.exp<-runSALSA2D(hfGLM,
                              salsa2dlist, 
                              d2k=distMatshavers$dataDist,
                              k2k=distMatshavers$knotDist,
                              basis = "exponential", ##
                              suppress.printout = FALSE)


# make predictions on response scale
nlPredData_havers$predsHavers.exp<-predict(newdata = nlPredData_havers,
                               g2k = preddist_havers,
                               object = salsa2dOutput$bestModel)

ggplot(nlPredData_havers[nlPredData$mindist<2500,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=(predsHavers.exp))) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()



################################################################
# Use variogram to get the central tendency
################################################################

#https://github.com/lindesaysh/MRSea/blob/HEAD/R/getRadiiSequence.R



rs_havers<-getRadiiSequence(method = "variogram",
                     numberofradii = 8,
                     xydata = hfNL[,c("x.pos", "y.pos")],
                     response = hfNL$response,
                     basis = "gaussian",
                     distMatrix = distMatshavers$dataDist)

rs.exp_havers<-getRadiiSequence(method = "variogram",
                         numberofradii = 8,
                         xydata = hfNL[,c("x.pos", "y.pos")],
                         response = hfNL$response,
                         basis = "exponential",
                         distMatrix = distMatshavers$dataDist)

# In this case, the range parameter is 10.09 which suggests on average, 
# the spatial correlation decays after approximately 4.5km. Which is equivalent
# to 18*log10(4567) or 65 dB




ggplot(nlPredData_havers[nlPredData$mindist<2369,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=(predsHavers.exp))) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()


ggplot(nlPredData_havers[nlPredData$mindist<2369,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=(preds))) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()



# Fit of variogram, gaussian
salsa2dlist$r_seq <- rs_havers
salsa2dOutput.vario.haversine <-runSALSA2D(model = hfGLM,
                                 salsa2dlist = salsa2dlist,
                                 d2k=distMatshavers$dataDist,
                                 k2k=distMatshavers$knotDist,
                                 basis = 'gaussian',
                                 suppress.printout = FALSE)


# make predictions on response scale
nlPredData_havers$predsHavers.vario<-predict(newdata = nlPredData_havers,
                                           g2k = preddist_havers,
                                           object = salsa2dOutput.vario.haversine$bestModel)


ggplot(nlPredData_havers[nlPredData$mindist<2369,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=(predsHavers.vario))) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()

# Fit of variogram, exponential
salsa2dlist$r_seq <- rs.exp_havers 
salsa2dOutput.vario.exp.haversine <-runSALSA2D(model = hfGLM,
                                           salsa2dlist = salsa2dlist,
                                           d2k=distMatshavers$dataDist,
                                           k2k=distMatshavers$knotDist,
                                           basis = 'exponential',
                                           suppress.printout = FALSE)


# make predictions on response scale
nlPredData_havers$predsHavers.vario.exp<-predict(newdata = nlPredData_havers,
                                             g2k = preddist_havers,
                                             object = salsa2dOutput.vario.exp.haversine$bestModel)


ggplot(nlPredData_havers[nlPredData$mindist<2369,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=(predsHavers.vario.exp))) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()