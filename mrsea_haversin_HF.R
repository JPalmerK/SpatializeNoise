# MRSEA modelling
rm(list = ls())
library(tidyr)
library(dplyr)
library(MRSea)
library(sp)
library(ggplot2)
library(lubridate)
library(geosphere)
library(ggcorrplot)
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
#csv_directory='F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023 Noise Files'
csv_directory='C:/Users/kaitlin.palmer/Documents/GitHub/NL_spatialize/SpatializeNoise/GPS_CSV-20230923T045356Z-001//MorroBay Mar 2023 Noise Files/'

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

nlLong =noiseDf[,c(1,8:25)] %>% 
  pivot_longer(cols = TOL_500: TOL_20000,
               names_to = 'Frequency',names_prefix = 'TOL_',
               values_to = 'metric')
nlLong$Dummy = as.factor(nlLong$Frequency)

ggplot(nlLong)+
  facet_wrap(~DriftName,ncol = 1)+
  geom_tile(aes(x=UTC, y= Frequency, fill = metric))

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

meanTimes = as.numeric(mean(allNl$UTC))

# Detrend the data with rolling median
allNl <- allNl %>%
  group_by(Band) %>%
  dplyr::arrange(desc(UTC)) %>% 
  dplyr::mutate(NL_ave = zoo::rollmedian(NL, k = 211, fill = NA)) %>%
  dplyr::mutate(Time =  as.numeric(UTC-meanTimes)) %>%
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


# Plot data and trends       
ggplot(allNl)+
  geom_line(aes(x=UTC, y=response, color = DriftName))+
  scale_color_brewer(palette = 'Paired', name= 'Drift ID')+
  facet_wrap(~Band, nrow = 2)+
  theme_bw()



# restrict the analysis to periods when all 8 buoys were in the water
dataStart = noise500reshape$UTC[100]
dataStop = noise500reshape$UTC[3561]

allNl= subset(allNl, UTC>=dataStart & UTC<= dataStop & Lat >35.3)

lfNL = subset(allNl, Band == '500Hz')
hfNL = subset(allNl, Band == '20kHz')


########################################################################
# Explore temporal autocorrelation
########################################################################

noiseReshape <- subset(allNl, Band == '20kHz') %>%
  spread(key = DriftName, value = NL)

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



rm(list=setdiff(ls(), c("allNl", 'lfNL', 'hfNL', 'GPSdf')))


################################################################
# Set up the knog grid, distance matrix, and prediction matrix
################################################################
source("~/GitHub/NL_spatialize/SpatializeNoise/makeDistsHaversine.R")
source("~/GitHub/NL_spatialize/SpatializeNoise/getRadiiSequence.R")

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

# Prediction Space
nlPredData_havers = expand.grid(Lon = 
                                  seq(min(allNl$Lon), max(allNl$Lon), 
                                      length.out =150),
                                Lat = 
                                  seq(min(allNl$Lat), max(allNl$Lat), 
                                      length.out =300))
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

########################################################################
# Fit the MRSEA model
########################################################################
# knotgridhavers_converted =as.matrix(coordinates(cord.UTM))
# colnames(knotgridhavers_converted)

medTime = median(hfNL$UTC)
hfNL$Time = as.numeric(hfNL$UTC- medTime)

# hfGLM <- glm(NL ~ 1, data=hfNL)
# hfGLM <- glm(NL ~ 1+ offset(NL_ave), data = hfNL)#***
# hfGLM <- glm(NL ~ 1+ offset(Time) , data = hfNL)#*
# hfGLM <- glm(response ~ 1+ offset(Time) , data = hfNL)#*
# hfGLM <- glm(response ~ 1+ offset(NL_ave) , data = hfNL)
# hfGLM <- glm(NL ~ 1+ offset(log(NL_ave)) , data = hfNL)
# hfGLM <- glm(NL ~ 1+ offset(log10(NL_ave)) , data = hfNL)
# hfGLM <- glm(response ~ 1+ offset(log10(NL_ave)), data=hfNL)#** 
# hfGLM <- glm(response ~ 1+ offset(log(NL_ave)), data=hfNL)#xxxx
hfGLM <- glm(response ~ 1, data=hfNL)#***



salsa2dlist<-list(fitnessMeasure = 'BIC',
                  knotgrid = knotgridhavers,
                  startKnots=20,
                  minKnots=10,
                  maxKnots=40,
                  splines = c("ns"))
#minKnots_1d = rep(1, length(varlist)))


# Gaussian basis
salsa2dOutput.haversine<-runSALSA2D(hfGLM,
                        salsa2dlist, 
                          d2k=distMatshavers$dataDist,
                          k2k=distMatshavers$knotDist,
                          basis = "gaussian", ##
                          suppress.printout = FALSE)

summary(salsa2dOutput.haversine$bestModel)
runDiagnostics(salsa2dOutput.haversine$bestModel)

cv.gamMRSea(hfNL, salsa2dOutput.haversine$bestModel, K=10)$delta[2]

nlPredData_havers$Time = median(hfNL$Time)
nlPredData_havers$NL_ave =median(hfNL$Median_NL)

# make predictions on response scale
preds<-predict(newdata = nlPredData_havers,
               g2k = preddist_havers,
               object = salsa2dOutput.haversine$bestModel)


nlPredData_havers$preds <-preds[,1]
nlPredData_havers$mindist =(apply(X = preddist_havers, 1, min))


animalsStart = as.POSIXct('2023-03-15 00:00:00', tz= 'UTC')
animalsStop = as.POSIXct('2023-03-15 10:00:00', tz= 'UTC')

ggplot(nlPredData_havers[nlPredData_havers$mindist<5000,]) +
  geom_tile(aes(x=Lon, y=Lat, fill=preds)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Easting (km)") + ylab("Northing (km)") + theme_bw()+
  ggtitle('High Frequency (20khz) Noise')+
  # geom_path(data =subset(GPSdf, UTC> animalsStart & UTC<animalsStop),
  #            aes(x=Longitude, y=Latitude, group= DriftName),
  #           color ='black')
  geom_point(data = hfNL, aes(x= Lon, y=Lat, color = response))




################################################################
# Use variogram to get the central tendency
################################################################

#https://github.com/lindesaysh/MRSea/blob/HEAD/R/getRadiiSequence.R



rs_havers<-getRadiiSequence(method = "variogram",
                     numberofradii = 15,
                     xydata = hfNL[,c("x.pos", "y.pos")],
                     response = hfNL$response,
                     basis = "gaussian",
                     distMatrix = distMatshavers$dataDist)

# In this case, the range parameter is 10.09 which suggests on average, 
# the spatial correlation decays after approximately 4.5km. Which is equivalent
# to 18*log10(4567) or 65 dB


ggplot(nlPredData_havers[nlPredData_havers$mindist<attr(rs_havers, 'vg.fit')[[3]][[2]],]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=preds)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation",
                       limits=c(-4, 6)) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_path(data = lfNL, aes(y=Lat, x=Lon, group = DriftName), 
            color = 'gray40')+
  ggtitle('Gaussian High Frequency (500Hz)')+
  theme_bw()


figData1 = allNl[,c('Lat', 'Lon', 'Band','NL')]
figData2 = allNl[,c('Lat', 'Lon', 'Band','response')]
colnames(figData1)[4]<-'Noise Metric'
colnames(figData2)[4]<-'Noise Metric'

figData = rbind(figData1, figData2) 
figData$panel = paste0(figData$Band, figData$`Noise Metric`)


ggplot(figData)+
  geom_point(aes(y=Lat, x=Lon, color = 'Noise Metric'))+
  facet_wrap(~panel, nrow =2)+
  scale_color_distiller(palette = "Spectral",name="Noise Level (dB)")+
  theme_bw()


# Fit of variogram, gaussian
salsa2dlist$r_seq <- rs_havers
salsa2dOutput.vario.haversine <-runSALSA2D(model = hfGLM,
                                 salsa2dlist = salsa2dlist,
                                 d2k=distMatshavers$dataDist,
                                 k2k=distMatshavers$knotDist,
                                 basis = 'gaussian',
                                 suppress.printout = FALSE)


runDiagnostics(salsa2dOutput.vario.haversine$bestModel)

# make predictions on response scale
nlPredData_havers$predsHavers.vario<-predict(newdata = nlPredData_havers,
                                           g2k = preddist_havers,
                                           object = salsa2dOutput.vario.haversine$bestModel)


ggplot(nlPredData_havers[nlPredData_havers$mindist<2760.64,]) + 
  geom_tile(aes(x=Lon, y=Lat, fill=predsHavers.vario)) +
  scale_fill_distiller(palette = "Spectral",name="NL Variation") +
  xlab("Longitude") + ylab("Latitude") + theme_bw()+
  ggtitle('Gaussian: variogram fit')



#########################################################################
# Cross validation
###########################################################################

# Kfold cross validation
xval = cv.gamMRSea(
  data =hfNL,
  modelobject =salsa2dOutput.haversine$bestModel,
  cost = function(y, yhat) mean((y - yhat)^2),
  K = 5,
  replicate = FALSE)

# CV scores
cv.gamMRSea(hfNL, salsa2dOutput.haversine$bestModel, K=5)$delta[2]
cv.gamMRSea(hfNL, salsa2dOutput.vario.haversine$bestModel, K=5)$delta[2]
# Variogram method has a lower CV