# Make plots of area monitored
rm(list =ls())


######################################################################
# Load and prep the data
#################################################################
library(sp)
library(ncdf4)
library(ggplot2)
library(lubridate)
library(dplyr)

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


# Add noise level data
csv_directory='F:\\GPS_CSV-20230923T045356Z-001\\MorroBay Mar 2023 Noise Files'
csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
dataframes_list <- list()

# Loop through each CSV file load, change name
for (csv_file in csv_files) {
  # Read the CSV file
  df <- read.csv(csv_file, header = TRUE)
  
  # Extract the first and eighth columns
  df <- df[, c(1, 2, 8,11,  21, 24)]
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
noiseDf$datetime_posix <- as.POSIXct(noiseDf$UTC, 
                                     format = "%Y-%m-%dT%H:%M:%OSZ",
                                     tz='UTC')

# Clean out data for drifts that don't have gps or noise levels
GPSdf= subset(GPSdf, DriftName %in% noiseDf$DriftName)
noiseDf= subset(noiseDf, DriftName %in% GPSdf$DriftName)

# Load the wind lease area
library(sf)
WLA <- st_read(
  "F:\\GPS_CSV-20230923T045356Z-001\\MorroBay_WEA_2021_11_12.shp")
WLA = as.data.frame(WLA[[6]][[1]][1])
colnames(WLA)<-c('UTMx','UTMy')

# Add lat/lon
# Add UTC coords with the sp package
cord.dec = SpatialPoints(cbind(WLA$UTMx,WLA$UTMy),
                         proj4string=CRS("+init=epsg:32610"))


cord.dec = spTransform(cord.dec, CRS("+proj=longlat"))
WLA$Lon =  coordinates(cord.dec)[,1]
WLA$Lat =  coordinates(cord.dec)[,2]

# Load the bthymetry
# Use the extents to download the GEBCo data from here:https://download.gebco.net/ 
nc_data <- nc_open('F:\\GPS_CSV-20230923T045356Z-001\\gebco_2023_n36.0791_s34.6025_w-122.8052_e-120.8013.nc')
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
elevation <- ncvar_get(nc_data, "elevation")
bathyRelief = data.frame(Lat=sort(x=rep(lat, length(lon))),
                         Long=rep(lon, length(lat)),
                         depth=c((elevation)))

# Add UTC coords with the sp package
cord.dec = SpatialPoints(cbind(bathyRelief$Long,bathyRelief$Lat),
                         proj4string=CRS("+proj=longlat"))
cord.dec = spTransform(cord.dec, CRS("+init=epsg:32610"))
bathyRelief$UTMx =  coordinates(cord.dec)[,1]
bathyRelief$UTMy =  coordinates(cord.dec)[,2]



# Check
p<-ggplot()+
  geom_path(data = GPSdf, aes(x=UTMx, y=UTMy, group = DriftName), color='gray30')+
  geom_path(data= WLA, aes(x= UTMx, y=UTMy), color = 'gold')+
  theme_bw()
p



############################################################################
# TL functions
############################################################################


# frequency dependent acoustic absorption
AcousticAbsorption <-function(f, Z=0, Temp=5, S=35, pH=8){
  
  # Following Kinsler, et al "Fundamentals of Acoustics, Fourth Edition" p. 226-228.
  # f = frequency in Hz
  # Z = depth in km
  # Temp = temperature in C
  # S = salinity in ppt
  # pH = pH
  
  
  f_1 = 780*exp(Temp/29)
  f_2 = 42000*exp(Temp/18)
  A = 0.083*(S/35)*exp(Temp/31 - Z/91 + 1.8*(pH-8))
  B = 22*(S/35)*exp(Temp/14-Z/6)
  C = 4.9e-10*exp(-Temp/26 - Z/25)
  boric_acid = A/(f_1^2+f^2) # contribution from boric acid
  MgSO4 = B/(f_2^2+f^2) # contribution from MgSO4
  hydrostatic = C # contribution from hydrostatic pressure
  alpha = (boric_acid + MgSO4 + hydrostatic)*f^2 #db/KM
  
  return(mean(alpha))
}

# Sonar equation over which to optimize r given, SL, NL, depth, frequency
# and SNR threshold
logfun <- function(SL, NL, SNRthresh, h, f) {
  function(r) {
    alpha = AcousticAbsorption(f)
    abs(SL - (15*log10(r)+10*log10(h/2)+(alpha/1000)*r) - NL - SNRthresh)
  }
}


# source, level, noise level, frequen(ies), and SNR threshold
f=500:900
SL_max=188
SL_min = 165
NL=80
SNRthresh=2
h=20


# Transmission loss
alpha = AcousticAbsorption(200)
TL<- function(r, h=100, alpha=0){
  (15*log10(r)+10*log10(h/2)+(alpha/1000)*r)}


#########################################################################
# Haversine fx
haversine_dist <- function(long1, lat1, long2, lat2) {
  rad <- pi/180
  a1 <- lat1*rad
  a2 <- long1*rad
  b1 <- lat2*rad
  b2 <- long2*rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1)*cos(b1)*(sin(dlon/2))^2
  c <- 2*atan2(sqrt(a), sqrt(1 - a))
  R <- 6378137
  d <- R*c
  return(d)
}
###############################################
# Create the time and noise level grid
###############################################

predGrid = expand.grid(DriftName =  unique(GPSdf$DriftName),
                       UTC = seq(min(GPSdf$UTC), max(GPSdf$UTC), by= '2 min'))
predGrid$Lat=NaN
predGrid$Lon =NaN
predGrid$NL = NaN
#################################################

# For each timestep, each drift, and each noise level band, get the area area monitored

measurements = colnames(noiseDf)[2:7]
drifters = unique(GPSdf$DriftName)
centerFs = c(125,500,1000,10000,20000)


dataOut = data.frame()


for(m in 1:length(measurements)){
  
  NLsub = noiseDf[,c(1, m+1, 8)]
  predGrid = expand.grid(UTC = seq(min(GPSdf$UTC),
                                   max(GPSdf$UTC), by= '2 min'))
  predGrid$Lat=NaN
  predGrid$Lon =NaN
  predGrid$NL = NaN
  predGrid$Band = measurements[m]

  
  
  alpha = AcousticAbsorption(centerFs[m])
  
  
  for(drift in drifters){
    
    predGrid$DriftName = drift
    
    GPSsub = subset(GPSdf, DriftName == drift)
    NLsub =  subset(noiseDf, DriftName==drift)
    
    UTMflon <- approxfun(GPSsub$UTC, GPSsub$Longitude)
    UTMflat <- approxfun(GPSsub$UTC, GPSsub$Latitude)
    NL_FUN <-  approxfun(NLsub$datetime_posix, NLsub[,2])
    
    predGrid$NL= NL_FUN(predGrid$UTC)
    predGrid$Lat = UTMflat(predGrid$UTC)
    predGrid$Lon = UTMflon(predGrid$UTC)
    
    print(drift)
    
    tStamps = which(!is.na(predGrid$NL))
  
    for(t in tStamps){
      
      # Preallocate the grid
      # Create 8 grids to show where the whale might be
      whaleGrid = expand.grid(
        UTMx= seq(min(GPSdf$UTMx)-2000, max(WLA$UTMx)+2000, by=1000),
        UTMy= seq(min(GPSdf$UTMy)-1000, max(WLA$UTMy)+3000, by=1000))
        DriftName = drift
      
      # Add lat/lon
      # Add UTC coords with the sp package
      cord.dec = SpatialPoints(cbind(whaleGrid$UTMx,whaleGrid$UTMy),
                               proj4string=CRS("+init=epsg:32610"))
      
      
      cord.dec = spTransform(cord.dec, CRS("+proj=longlat"))
      whaleGrid$Lon =  coordinates(cord.dec)[,1]
      whaleGrid$Lat =  coordinates(cord.dec)[,2]
      whaleGrid$cellId =paste0(whaleGrid$UTMx, whaleGrid$UTMy)
      whaleGrid$NL = NL_FUN(predGrid$UTC[t])
      
      # Range from each point to each Drift
      # Likely an overestimate of the error
      whaleGrid$rElipsoid = 
        haversine_dist(whaleGrid$Lon,
                       whaleGrid$Lat,
                       predGrid$Lon[t], 
                       predGrid$Lat[t])
      
      # The expected SNR grids
      whaleGrid$ExpectedSNR_havers = 177-whaleGrid$NL-TL(whaleGrid$rElipsoid, alpha)
      
    }
    
  }
  
}




