########################################################
# Simulation of SNR and TDOA
################################################
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
  df <- df[, c(1, 8)]
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



# Check
p<-ggplot()+
  geom_path(data = GPSdf, aes(x=UTMx, y=UTMy, group = DriftName), color='gray30')+
  geom_path(data= WLA, aes(x= UTMx, y=UTMy), color = 'gold')+
  theme_bw()
p


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



rm(list=setdiff(ls(), c("GPSdf", "WLA","noiseDf", 'p','bathyRelief')))


#########################################################################
# Create a moving animal in the space
#########################################################################

library(SiMRiv)

#set.seed(20)

sim.crw = data.frame()

# Make a bunch of whale simulations
for(ii in 1:20){
  # Simulate call source levels and times
  MinTime =  as.POSIXct(min(GPSdf$UTCDatehour), tz = 'UTC')+hours(12)
  MaxTime =   as.POSIXct(max(GPSdf$UTCDatehour), tz = 'UTC')-hours(24)
  
  # Create the random movement and center it
  c.rand.walker <- species(state.CRW(.99)+15)
  sim.set <- as.data.frame(simulate(c.rand.walker, 800))
  
  sim.set$UTMx=sim.set$V2+runif(1, min(GPSdf$UTMx), max(WLA$UTMx))
  sim.set$UTMy=sim.set$V1+runif(1,min(WLA$UTMy-20000) ,max(WLA$UTMy))
  sim.set$UTC = seq(MinTime, MaxTime, length.out = nrow(sim.set))
  
  sim.set$WhaleId= ii
  
  # simulate calling
  sim.set$SL =rbinom(n = nrow(sim.set),size = 1, prob = .2)*
    rnorm(n = nrow(sim.set), mean = 177, sd = 4)
  sim.set$SL[sim.set$SL==0]=NaN
  
  sim.crw=rbind(sim.crw, sim.set)
  
}

# Show the GPS tracks, whale, and exclusion area
#p = p+geom_path(data= sim.crw, aes(x=UTMx, y=UTMy, group=whaleGrid), color='red')



# check the plot
p= p+geom_point(data= subset(sim.crw, !is.na(SL)), aes(x=UTMx, y= UTMy, 
                                                       size = SL, color=UTC))

p
# 
# p +geom_path(data = subset(GPSdf, UTC >= min(sim.crw$UTC) & UTC <= max(sim.crw$UTC)),
#              aes(x = UTMx, y=UTMy, group= DriftName), color= 'cadetblue4',
#              size =2)

# Add lat/lon
# Add UTC coords with the sp package
cord.dec = SpatialPoints(cbind(sim.crw$UTMx,sim.crw$UTMy),
                         proj4string=CRS("+init=epsg:32610"))


cord.dec = spTransform(cord.dec, CRS("+proj=longlat"))
sim.crw$Lon =  coordinates(cord.dec)[,1]
sim.crw$Lat =  coordinates(cord.dec)[,2]


rm(list=setdiff(ls(), c("GPSdf", "WLA","noiseDf", 'p', 'sim.crw',
                        'alpha', 'TL','bathyRelief')))

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
alpha = AcousticAbsorption(800)
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

###########################################################################
# Calculate the recieved levels on each hydrophone


# clear out the data without source levels 
sim.calls=sim.crw[!is.na(sim.crw$SL),]
sim.calls$inWLA = point.in.polygon(sim.calls$UTMx, sim.calls$UTMy,
                                   WLA$UTMx, WLA$UTMy)
sim.calls$propAreainWLA = NaN
rownames(sim.calls)<-1:nrow(sim.calls)

PositionError = 1152/1500 # 80th percentile of the GPS error
SSPError = (1500*.2)/1500 # 20% error in speed of sound
TDOA_error = SSPError+PositionError

figs =list()
sim.calls$Blobn
sim.calls$Blob
sim.calls$Correct =FALSE

colnames(bathyRelief)[c(1,2)]<-c('Latitude', 'Longitude')
#Basemap
p<-ggplot()+
  geom_tile(data= bathyRelief, aes(x = Longitude, y=Latitude, fill= depth))+
  geom_path(data= WLA, aes(x= Lon, y=Lat), color = 'gold',linewidth=2)+
  xlim(c(-122.5, -121.1))+
  ylim(c(35.25,36))+
  theme_bw()


sim.calls$MaxSNR=NaN
sim.calls$GridCenterLat=NaN
sim.calls$GridCenterLon=NaN
sim.calls$rNearestInstrument=NaN
sim.calls$NearestInstrumentLat =NaN
sim.calls$NearestInstrumentLon =NaN
sim.calls$AllDetInstLat= NaN
sim.calls$AllDetInstLon= NaN
sim.calls$AllDetInstID = NaN

acceptedLocs = list()
whaleGridTemplate = expand.grid(
  UTMx= seq(min(WLA$UTMx)-10000, max(WLA$UTMx)+10000, by=1000),
  UTMy= seq(min(GPSdf$UTMy)-1000, max(WLA$UTMy)+10000, by=1000),
  DriftName = unique(GPSdf$DriftName))


# Add lat/lon
# Add UTC coords with the sp package
cord.dec = SpatialPoints(cbind(whaleGridTemplate$UTMx,whaleGridTemplate$UTMy),
                         proj4string=CRS("+init=epsg:32610"))


cord.dec = spTransform(cord.dec, CRS("+proj=longlat"))
whaleGridTemplate$Lon =  coordinates(cord.dec)[,1]
whaleGridTemplate$Lat =  coordinates(cord.dec)[,2]
whaleGridTemplate$cellId =paste0(whaleGridTemplate$UTMx, whaleGridTemplate$UTMy)

# Range from each point to each Drift
whaleGridTemplate$rElipsoid=NaN
whaleGridTemplate$SNRok_havers=NaN

DriftNames = unique(noiseDf$DriftName)

# Step through each call, determine acceptable SNR and TDOA ranges
#for(ii in 1:nrow(sim.calls)){
for(ii in 1:5){ 
  
  # Preallocate the grid
  # Create 8 grids to show where the whale might be
  whaleGrid = whaleGridTemplate
  
  
  # Recieved level dataframe
  RLdf = data.frame(DriftName = as.factor(DriftNames))
  
  for(drift in DriftNames){
    GPSsub = subset(GPSdf, DriftName == drift)
    NLsub =  subset(noiseDf, DriftName==drift)
    
    UTMflon <- approxfun(GPSsub$UTC, GPSsub$Longitude)
    UTMflat <- approxfun(GPSsub$UTC, GPSsub$Latitude)
    
    NL_FUN <-  approxfun(NLsub$datetime_posix, NLsub$TOL_500)
    
    ##############################################################
    # Calculate the range from the whale to the GPS, TDOA and RL
    ############################################################
    
    # Lat/lon/ of the drift when the call was produced
    RLdf$Lon[RLdf$DriftName==drift] =UTMflon(sim.calls$UTC[ii])
    RLdf$Lat[RLdf$DriftName==drift] = UTMflat(sim.calls$UTC[ii])
    RLdf$NoiseLevel[RLdf$DriftName==drift] = NL_FUN(sim.calls$UTC[ii])
    
    
    # Spatial Uncertainty, GPS- estimate the speed at a the call arrival time
    # then how far it could drift
    idxMindiff = which.min(abs((sim.calls$UTC[ii])- (GPSsub$UTC)))
    Time_drift = seconds(GPSsub$UTC[idxMindiff]-seconds(sim.calls$UTC[ii]))
    
    # Likely an overestimate of the error
    whaleGrid$rElipsoid[whaleGrid$DriftName==drift] = 
      haversine_dist(whaleGrid$Lon[whaleGrid$DriftName==drift],
                     whaleGrid$Lat[whaleGrid$DriftName==drift],
                     RLdf$Lon[RLdf$DriftName==drift], 
                     RLdf$Lat[RLdf$DriftName==drift])
  }
  
  # Range from the whale to the hydrophone in kms  (truth)
  RLdf$range_havers = haversine_dist(sim.calls$Lon[ii],
                                     sim.calls$Lat[ii],
                                     RLdf$Lon, RLdf$Lat)
  
  RLdf$ArrialTimeHavers =  sim.calls$UTC[ii]+(RLdf$range_havers/1500)
  RLdf$SNRHavers= sim.calls$SL[ii]-TL(RLdf$range_havers)-RLdf$NoiseLevel
  
  # Call detected or not
  RLdf$Detected = ifelse(RLdf$SNRHavers>12,1,0)
  RLdf$DetecteColor = ifelse(RLdf$SNRHavers>12,'green','maroon')
  
  
  # estimate the drifter location from the available data
  whaleGrid= merge(whaleGrid, 
                   RLdf[,c('ArrialTimeHavers','DriftName',
                           'SNRHavers', 'NoiseLevel', 'Detected')], 
                   by.y='DriftName')
  
  # The expected SNR at each grid location- mean
  whaleGrid$ExpectedSNR_havers = 177-TL(whaleGrid$rElipsoid)-
    whaleGrid$NoiseLevel
  
  # For detected calls
  whaleGrid$SNRok_havers[whaleGrid$Detected==1]=
    ifelse(whaleGrid$ExpectedSNR_havers[whaleGrid$Detected==1]+11  >= 
             whaleGrid$SNRHavers[whaleGrid$Detected==1] &
             (whaleGrid$ExpectedSNR_havers[whaleGrid$Detected==1]-12) <= 
             whaleGrid$SNRHavers[whaleGrid$Detected==1], 1,0)
  

  # 
   # For not detected calls the SNR must be lower than the minimum expected SNR
   whaleGrid$SNRok_havers[whaleGrid$Detected==0]=
     ifelse(whaleGrid$ExpectedSNR_havers[whaleGrid$Detected==0]<=24,1,0)
    
  # # which points matched the arrival SNR
  # whaleGrid$SNRok_havers=
  #   ifelse(whaleGrid$ExpectedSNR_havers+11  >= whaleGrid$SNRHavers &
  #            (whaleGrid$ExpectedSNR_havers-12) <= whaleGrid$SNRHavers, 1,0)
  # 
  # # For non-detected calls, keep all SNR below 12 db
  # idxWhale = which(whaleGrid$Detected==0)
  # whaleGrid$SNRok_havers[idxWhale]= ifelse(whaleGrid$ExpectedSNR_havers[idxWhale]>12,0,1)
  
  
  # Count up the cells that were OK by all buoys
  aggData = aggregate(whaleGrid, SNRok_havers~cellId,FUN=sum)
  colnames(aggData)[2]<-'Count'
  
  # Which grid locs were 'ok' by all receivers?
  whaleGrid = merge(whaleGrid, aggData, by = 'cellId', all.x=TRUE)
  
  #RLDFDetected = subset(RLdf, Detected==1)
  
  # Pull out the accepted locations, usefule for all calls even if for units
  # where animal wasn't detected also nix anything more than 50kms out
  accepted_locs = subset(whaleGrid, Count == nrow(RLdf))
  
  ########################################################
  # Create n-1 tdoa trids, nix data outside the expected TDOA values
  ########################################################
  if(sum(RLdf$Detected)>1){
    
    # All TDOA combinations
    combinations = combn(RLdf$DriftName[RLdf$Detected==1], 2,simplify = TRUE)
    
    for(kk in 1:ncol(combinations)){
      
      # Pull out the grids for each drift
      grid1 = subset(accepted_locs, DriftName==combinations[1, kk])
      grid2 = subset(accepted_locs, DriftName==combinations[2, kk])
      
      # Caluclate the expected TDOA
      grid1$TDOA = (grid1$rElipsoid-grid2$rElipsoid)/1500
      
      
      #Observed TDOA between the pairs
      OBS_tdoa = RLdf$ArrialTimeHavers[RLdf$DriftName==combinations[1,kk]]-
        RLdf$ArrialTimeHavers[RLdf$DriftName==combinations[2,kk]]
      
      
      # Pick out the locations that are OK within the TDOA grid
      grid1$ok =(grid1$TDOA>(OBS_tdoa-TDOA_error) & 
                   grid1$TDOA<=(OBS_tdoa+TDOA_error))
      
      # Trim the accepted locations
      accepted_locs= subset(accepted_locs, 
                            cellId %in% grid1$cellId[grid1$ok==TRUE])
      
    }
    # Add the SNR, and grid center
    sim.calls$MaxSNR[ii]=max(RLdf$SNRHavers[RLdf$Detected==1])
    sim.calls$GridCenterLat[ii]=min(accepted_locs$Lat)+
      diff(abs(range(accepted_locs$Lat)))/2
    sim.calls$GridCenterLon[ii]=min(accepted_locs$Lon)+
      diff(abs(range(accepted_locs$Lon)))/2
    
    # Distance to the nearest instrument
    sim.calls$rNearestInstrument[ii]=min(RLdf$range_havers[RLdf$Detected==1])
    sim.calls$NearestInstrumentLat[ii] = RLdf$Lat[which.max(RLdf$SNRHavers)]
    sim.calls$NearestInstrumentLon[ii] = RLdf$Lon[which.max(RLdf$SNRHavers)]
    
    # Distance to all instruments
    sim.calls$AllDetInstLat[ii] =  list(RLdf$Lat[RLdf$Detected==1])
    sim.calls$AllDetInstLon[ii] = list(RLdf$Lon[RLdf$Detected==1])
    sim.calls$AllDetInstID[ii] = list(RLdf$DriftName[RLdf$Detected==1])
    
    # Create a figure
    figs[[ii]]=p+
      geom_point(data = accepted_locs, aes(x=Lon, y=Lat), color = 'gray60')+
      geom_path(data = GPSdf[GPSdf$UTC<= sim.calls$UTC[ii],], 
                aes(x=Longitude, y=Latitude, group = DriftName), color='black')+
      geom_point(data = RLdf, 
                 aes(x=Lon, y=Lat, group = DriftName, shape = DriftName, 
                     color= as.factor(Detected)))+
      geom_path(data= sim.crw[sim.crw$UTC<=sim.calls$UTC[ii],], 
                aes(x=Lon, y=Lat, group=WhaleId), color='red')+
      geom_point(data=sim.calls[ii,], 
                 aes(Lon, Lat), size=3, color='red')
    
    # Save the accepted locations
    acceptedLocs[[ii]]<-accepted_locs
  }
  
  
  print(ii)
}

## Clean out non-detected calls
sim.calls= sim.calls[!is.na(sim.calls$GridCenterLat),]

# Get the distance between the grid center and the whale
sim.calls$rgrid_to_whale = haversine_dist(sim.calls$Lon, sim.calls$Lat,
                                          sim.calls$GridCenterLon,
                                          sim.calls$GridCenterLat)




# 
# p+
#   geom_path(data = GPSdf, 
#             aes(x=Longitude, y=Latitude, group = DriftName), color='black')+
#   geom_path(data= sim.crw, 
#             aes(x=Lon, y=Lat, group = WhaleId), color='red')+
#   geom_point(data=sim.calls, 
#              aes(Lonitude, Latitude), size=2)+
#   geom_point(data= sim.calls, 
#              aes(GridCenterLon, GridCenterLat), color='green')
#   xlim(c(-122.5,-121.25))+
#   ylim(c(35, 36))+ 
#   theme_bw()


  
  #############################################################
  # Min Critters
  ############################################################
  
  # ii =20
  # jj=240
  # 
  # aa = do.call(cbind.data.frame, acceptedLocs[[ii]])
  # bb = do.call(cbind.data.frame, acceptedLocs[[jj]])
  # 
  # whaleGrid$Overlap = ifelse(whaleGrid$cellId %in% aa$cellId, 1,0)
  # whaleGrid$Overlap = whaleGrid$Overlap+
  #   ifelse(whaleGrid$cellId %in% bb$cellId, 1,0)
  # newLocs = subset(whaleGrid, Overlap>0)
  # 
  # ii=2
  # p+
  #   geom_point(data = newLocs, aes(x=Lon, y=Lat, fill = Overlap))+
  #   geom_path(data = GPSdf[GPSdf$UTC<= sim.calls$UTC[jj],], 
  #             aes(x=Longitude, y=Latitude, group = DriftName), color='black')+
  #   geom_point(data=sim.calls[c(ii,jj),], 
  #              aes(Lon, Lat), size=3, color='red')
  
  
  
  
  

    
#############################################################
# Get ancillary data
############################################################

library(PAMpal)
colnames(sim.calls)[c(9,10)]<-c('Longitude', 'Latitude')

# Get the data at the instrument with the loudest SNR, the grid center,
# and the whale. 

sim.callsLoudestSNR = sim.calls[,c('NearestInstrumentLat',
                                   'NearestInstrumentLon', 'UTC', 'MaxSNR')]
sim.callsLoudestSNR$Treatment = 'InstLocLoudestSNR'
colnames(sim.callsLoudestSNR)[c(1,2)]<-c('Latitude','Longitude')
sim.callsGridCenter = sim.calls[,c('GridCenterLat', 'GridCenterLon',
                                   'UTC','MaxSNR')]
sim.callsGridCenter$Treatment = 'ProbGridCenter'
colnames(sim.callsGridCenter)[c(1,2)]<-c('Latitude','Longitude')

sim.callsWhale = sim.calls[,c('Latitude','Longitude', 'UTC', 'MaxSNR')]
sim.callsWhale$Treatment = 'Truth'


#Pull out the data from all the locations at the time of the whale
library(dplyr)
sim.calls$row_id = 1:nrow(sim.calls)

allBuoyLocs <- sim.calls %>%
  tidyr::unnest(AllDetInstLat, AllDetInstLon, AllDetInstID) %>%
  dplyr::select(AllDetInstLat, AllDetInstLon, AllDetInstID, row_id, UTC)

colnames(allBuoyLocs)<- c('Latitude', 'Longitude','DriftName', "row_id", 'UTC')
allBuoyLocs$Treatment = 'AllBuoys'
allBuoyLocs$MaxSNR = NaN

allBuoyLocs = rbind(allBuoyLocs, 
                    cbind(sim.callsWhale, 
                          data.frame(row_id = 1:nrow(sim.callsWhale),
                                     DriftName = rep('Whale', nrow(sim.callsWhale)))))

#########################################################################
# Get environmental Data
#########################################################################

# Evaluation data
evalData = rbind(sim.callsLoudestSNR, sim.callsGridCenter, sim.callsWhale)
evalData<-matchEnvData(evalData,nc='erdMBchla8day', var=c('chlorophyll'))
allBuoyLocs<-matchEnvData(allBuoyLocs,nc='erdMBchla8day', var=c('chlorophyll'))


evalData<-matchEnvData(evalData,nc='jplMURSST41mday', var=c('sst'))
allBuoyLocs<-matchEnvData(allBuoyLocs,nc='jplMURSST41mday', var=c('sst'))


evalData<-matchEnvData(evalData,nc='erdMH1pp8day', var=c('productivity'))
allBuoyLocs<-matchEnvData(allBuoyLocs,nc='erdMH1pp8day', var=c('productivity'))


evalData<-matchEnvData(evalData,nc='erdSrtm30plusSeafloorGradient', 
                          var=c('magnitude_gradient','sea_floor_depth' ))
allBuoyLocs<-matchEnvData(allBuoyLocs,nc='erdSrtm30plusSeafloorGradient', 
                       var=c('magnitude_gradient','sea_floor_depth' ))

evalData<-matchEnvData(evalData,nc = getEdinfo()[['HYCOM']], 
                          var=c('salinity', 'water_u', 'water_v'),
                       timeout=360)
allBuoyLocs<-matchEnvData(allBuoyLocs,nc = getEdinfo()[['HYCOM']], 
                       var=c('salinity', 'water_u', 'water_v'),
                       timeout=360)


# Distance to shore
dist2shore <- erddapToEdinfo(dataset='dist2coast_1deg',
                             baseurl='https://pae-paha.pacioos.hawaii.edu/erddap/')

# Then use this as the "nc" argument
evalData <- matchEnvData(evalData, nc=dist2shore)
allBuoyLocs <- matchEnvData(allBuoyLocs, nc=dist2shore)


# For each variable of interest, compare the error
EnvVars = c('chlorophyll_mean', 'sst_mean', 'productivity_mean', 'dist_mean',
            'magnitude_gradient_mean','sea_floor_depth_mean')

errorDF = data.frame()

for(env in EnvVars){
  
  dataSub = evalData[, c('Treatment', env, 'MaxSNR')]
  
  drifterError = 
    abs(dataSub[dataSub$Treatment=='Truth',2])-
    abs(dataSub[dataSub$Treatment=='InstLocLoudestSNR',2])
  
  blobError = 
    abs(dataSub[dataSub$Treatment=='Truth',2])-
    abs(dataSub[dataSub$Treatment=='ProbGridCenter',2])
  
  # Difference in error for each point
  errorDF = rbind(errorDF, 
                    data.frame(MaxSNR=dataSub$MaxSNR,
                      ErrorDiff = drifterError- blobError,
                               EnvVar = rep(env, length(blobError))))
  
}


# SNR bins 
errorDF$SNR_bins <- cut(errorDF$MaxSNR, 
                               breaks = seq(min(errorDF$MaxSNR),
                                            max(errorDF$MaxSNR), 
                                            length.out = 6), 
                               include.lowest = TRUE)

# Create the boxplot using ggplot2
ggplot(errorDF, aes(x = SNR_bins, y = ErrorDiff)) +
  facet_wrap(~EnvVar, scales = 'free_y')+
  geom_boxplot() +
  labs(x = "SNR Bins", y = "Difference in Measurement Error") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+ 
  theme_minimal()



  #############################################################
  # Getbackgorund levels
  ############################################################  
  
  colnames(bathyRelief)[1:2]
  colnames(bathyRelief)[1:2]<-c('Latitude', 'Longitude')
  bathyRelief$UTC<- median(GPSdf$UTC)
  
  bathyRelief<-matchEnvData(bathyRelief,nc='erdMBchla8day', var=c('chlorophyll'))
  bathyRelief<-matchEnvData(bathyRelief,nc='jplMURSST41mday', var=c('sst'))
  bathyRelief<-matchEnvData(bathyRelief,nc='erdMH1pp8day', var=c('productivity'))
  bathyRelief<-matchEnvData(bathyRelief,nc='erdSrtm30plusSeafloorGradient', 
                            var=c('magnitude_gradient','sea_floor_depth' ))
  bathyRelief<-matchEnvData(bathyRelief,evalData,nc = getEdinfo()[['HYCOM']], 
                            var=c('salinity', 'water_u', 'water_v'))
  
  # Distance to shore
  dist2shore <- erddapToEdinfo(dataset='dist2coast_1deg',
                               baseurl='https://pae-paha.pacioos.hawaii.edu/erddap/')
  
  # Then use this as the "nc" argument
  bathyRelief <- matchEnvData(bathyRelief, nc=dist2shore)
  
  envPlots = list()
  for(ii in 1:length(EnvVars)){
    
    dataSub = bathyRelief[, c('Longitude', 'Latitude', EnvVars[ii])]
    
    envPlots[[ii]]= 
      ggplot(dataSub)+
      geom_tile(aes(x=Longitude, y=Latitude, fill = dataSub[,3]))+
      geom_path(data= WLA, aes(x= Lon, y=Lat), color = 'gold',linewidth=2)+ 
      geom_path(data = GPSdf,
                aes(x=Longitude, y=Latitude, group = DriftName), color='red')+
      geom_path(data = sim.calls, aes(Longitude, Latitude, 
                                      color = as.factor(WhaleId)), linewidth=2)+
      xlim(c(-122.5, -121.1))+
      ylim(c(35.25,36))+
      theme_bw()
    
    rm(list=ls(dataSub))
    print(EnvVars[ii])
    
  }
  
  
  
#############################################################
# Turn into gif
############################################################
# 
# 
# # Use the extents to download the GEBCo data from here:https://download.gebco.net/ 
# nc_data <- nc_open('F:\\GPS_CSV-20230923T045356Z-001\\gebco_2023_n36.0791_s34.6025_w-122.8052_e-120.8013.nc')
# lon <- ncvar_get(nc_data, "lon")
# lat <- ncvar_get(nc_data, "lat", verbose = F)
# elevation <- ncvar_get(nc_data, "elevation")
# bathyRelief = data.frame(Lat=sort(x=rep(lat, length(lon))),
#                          Long=rep(lon, length(lat)),
#                          depth=c((elevation)))
# 
# # Add UTC coords with the sp package
# cord.dec = SpatialPoints(cbind(bathyRelief$Long,bathyRelief$Lat),
#                          proj4string=CRS("+proj=longlat"))
# cord.dec = spTransform(cord.dec, CRS("+init=epsg:32610"))
# bathyRelief$UTMx =  coordinates(cord.dec)[,1]
# bathyRelief$UTMy =  coordinates(cord.dec)[,2]
# 
# ggplot(bathyRelief)+geom_tile(aes(x = UTMx, y=UTMy, fill= depth))
# ggplot(bathyRelief)+geom_tile(aes(x = Long, y=Lat, fill= depth))
# 
# 
# library(gganimate)
# # Create an animation from the list of ggplots
# animation <- plot_list(figs)
# 
# # Specify the filename and other parameters
# anim_save("your_animation.gif", animation, renderer = gifski_renderer())
# 
# 
# 
# 
# #############################################################
