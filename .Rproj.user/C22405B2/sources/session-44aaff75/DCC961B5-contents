# Functions for spatial prediction

library(akima)
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


rescale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
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
  alpha=mean(alpha/1000) # dB/m
  return(alpha)
}

# Sonar equation over which to optimize r given, SL, NL, depth, frequency
# and SNR threshold
logfun <- function(SL, NL, SNRthresh, h, f) {
  function(r) {
    alpha = AcousticAbsorption(f)
    abs(SL - (15*log10(r)+10*log10(h/2)+(alpha/1000)*r) - NL - SNRthresh)
  }
}


# # source, level, noise level, frequen(ies), and SNR threshold
# f=500:900
# SL_max=188
# SL_min = 165
# NL=80
# SNRthresh=2
# h=20



TL<- function(r,f=10:30, c=1500, g= 0.016, d =30, H=90, S=3){
  # Return the transmission loss following Urick pg 153 and 
  #  Au and Hastings 107, eq 4.66 
  # f, frequency in hz
  # c, speed of sound in meters
  # g, unknown in meters
  # d, source depth in meters
  # H, layer depth in meters
  # S, seastate
  
  alpha = AcousticAbsorption(f)
  R = c/g
  r0 = sqrt((R*H)/2)*sqrt(H/(H-d))
  alphaLeakage = 2*S*((mean(f)/1000)/H)
  TL = 10*log10(r0)+10*log10(r)+(alpha+alphaLeakage)*r
  return(TL)
}

# 2d interpolation



#########################################################
# Simulate whale locations based on underlying density grid
#########################################################

createSimWhales = function(denGridSub, f, SLmean, SLstdev, h){
  # input: Density grid with normalized density 'NormDen' column
  
  # Make whales and figure out times
  denGridSub$Finwhales = rbinom(prob =denGridSub$NormDen, 
                                size = 1, n = nrow(denGridSub))
  
  
  sim.calls = denGridSub[denGridSub$Finwhales==1, c('lat', 'lon')]
  sim.calls= sim.calls[!is.na(sim.calls$lat),]
  colnames(sim.calls)<-c('Lat', 'Lon')
  sim.calls$Lat = jitter(sim.calls$Lat,amount = .06)
  sim.calls$Lon = jitter(sim.calls$Lon, amount = .1)
  
   
  sim.calls$SL = rnorm(SLmean, SLstdev, n=nrow(sim.calls)) 

  return(sim.calls)
  
}


#######################################################
# Figure out which lat/lons are acceptable based on SNR
######################################################

acceptedSNR_AreaMonitored<-function(SL=c(165,195), SNRthresh=15, NL, 
                                    senLat, senLon, whaleGrid, 
                                    SNR_observed =NaN){
  # Function to calcuate 1) at which lat/lons in whalegrid a call of SL and
  # noise level of NL could be detected. SNR must be above SNRthresh
  # Input:
  # SL - min and max SL of the calling animal
  # Returns:
  # dataOut - dataframe with the probabiliy of detection at each cell (Pdet),
  # and which cells fell in the acceptable range (SNR_oK)
  
  # Dist between the drift and the prediction grid cells
  dataOut = whaleGrid[,c('Lat', 'Lon')]
  rElipsoid =  haversine_dist(dataOut$Lon, dataOut$Lat,senLon, senLat)
  
  # If a SNR is given, then figure out which locations are good
  # Expected SNR at each grid cell if the call was in the center and the whale
  # was yelling
  SNR_min = SL[1]-TL(rElipsoid)-NL
  
  
  # Expected SNR at each grid cell if the call was in the center and the whale
  # was whispering
  SNR_max = SL[2]-TL(rElipsoid)-NL
  
  # what was the probability of detection (dirac fx, 0 or 1)
  dataOut$Pdet = ifelse(SNR_max>SNRthresh,1,0)
  
  if(!is.nan(SNR_observed)){
    dataOut$SNRok_havers =0
    #Which grid cells are ok by both the Min and Max SNR. I.e. the SNR must be 
    # between the min and max SNR
    dataOut$SNRok_havers[dataOut$Pdet==1]=
      ifelse(dataOut$SNR_max   >= SNR_observed &
               (dataOut$SNR_min) <= SNR_observed , 1,0)
    
    dataOut$SNR_ok = ifelse( SNR_observed >= dataOut$SNR_min & 
                              SNR_observed <= dataOut$SNR_max, 1,0)
  }
  
  return(dataOut)
  
}


calcWhaleGrid<-function(SL, SNRthresh, drifterLoc, gridArea){
  
  SL_diff = diff(SL)

  # Create the base grid for all units
  whaleGrid = expand.grid(
    Lat = unique(gridArea$Lat),
    Lon = unique(gridArea$Lon),
    DriftName = unique(drifterLoc$DriftName))
  
  whaleGrid$cellId =paste0(whaleGrid$Lat, whaleGrid$Lon)
  
  colnames(drifterLoc)[c(3,4,6)]<- c('LatDrifter', 'LonDrifter', 'rAnimal')
  
  # Merge the drifter information with the whale grid
  whaleGrid = merge(whaleGrid, drifterLoc, by= 'DriftName')
  
  # Distance between the drifter lcation and the grid cells
  whaleGrid$rElipsoid = haversine_dist(whaleGrid$Lon, whaleGrid$Lat,
                                       whaleGrid$LonDrifter, whaleGrid$LatDrifter)
   
  # The maximum expected SNR at each grid location- mean
  whaleGrid$ExpectedSNR_havers = mean(SL)-TL(whaleGrid$rElipsoid)-
    whaleGrid$NL
  whaleGrid$SNRok_havers=NaN
  
  # For detected calls
  whaleGrid$SNRok_havers[whaleGrid$Detected==1]=
    ifelse(whaleGrid$ExpectedSNR_havers[whaleGrid$Detected==1]+(SL_diff/2)  >= 
             whaleGrid$SNR[whaleGrid$Detected==1] &
             (whaleGrid$ExpectedSNR_havers[whaleGrid$Detected==1]-SL_diff/2) <= 
             whaleGrid$SNR[whaleGrid$Detected==1], 1,0)
  
  
  # For not detected calls the SNR must be lower than the minimum expected SNR
  whaleGrid$SNRok_havers[whaleGrid$Detected==0]=
    ifelse(whaleGrid$ExpectedSNR_havers[whaleGrid$Detected==0]<=SNRthresh,1,0)
  
  # Count up the cells that were OK by all buoys
  aggData = aggregate(whaleGrid, SNRok_havers~cellId,FUN=sum)
  colnames(aggData)[2]<-'Count'
  
  # Which grid locs were 'ok' by all receivers?
  whaleGrid = merge(whaleGrid, aggData, by = 'cellId', all.x=TRUE)
  whaleGrid$SNR_accepted = ifelse(whaleGrid$Count==nrow(drifterLoc),1,0)
  
  # Pull out the accepted locations, usefule for all calls even if for units
  # where animal wasn't detected also nix anything more than 50kms out
  accepted_locs = subset(whaleGrid, SNR_accepted == 1)
  
  
  ########################################################
  # Create n-1 tdoa trids, nix data outside the expected TDOA values
  ########################################################
  if(sum(drifterLoc$Detected)>1){
    
    # All TDOA combinations
    combinations = combn(drifterLoc$DriftName[drifterLoc$Detected==1],
                         2,simplify = TRUE)
    
    for(kk in 1:ncol(combinations)){
      
      # Pull out the grids for each drift
      grid1 = subset(accepted_locs, DriftName==combinations[1, kk])
      grid2 = subset(accepted_locs, DriftName==combinations[2, kk])
      
      # Caluclate the expected TDOA
      grid1$TDOA = (grid1$rElipsoid-grid2$rElipsoid)/1500
      
      
      #Observed TDOA between the pairs
      OBS_tdoa = drifterLoc$CallArrivalTime[drifterLoc$DriftName==combinations[1,kk]]-
        drifterLoc$CallArrivalTime[drifterLoc$DriftName==combinations[2,kk]]
      
      
      # Pick out the locations that are OK within the TDOA grid
      grid1$ok =(grid1$TDOA>(OBS_tdoa-TDOA_error) & 
                   grid1$TDOA<=(OBS_tdoa+TDOA_error))
      
      # Trim the accepted locations
      accepted_locs= subset(accepted_locs, 
                            cellId %in% grid1$cellId[grid1$ok==TRUE])
      
    }
    
  }
  
  # Add the accepted locations back into the whale grid
  whaleGrid$TDOA_and_SNRok = ifelse(whaleGrid$cellId %in% accepted_locs$cellId,
                                    1,0)
  
  return(whaleGrid)
}
