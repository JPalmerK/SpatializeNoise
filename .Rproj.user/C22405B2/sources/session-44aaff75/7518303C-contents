rm(list = ls())
library(mgcv)
library(dplyr)
library(lubridate)

load('Adrift_GPS_2023.rds')
load('Adrift_NL_2023_500Hz.rds')

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

# Link the gps and the noise levels
noiseDf$Lon = NaN
noiseDf$Lat = NaN

DriftNames = unique(GPSdf$DriftName)

for(drift in DriftNames){
  GPSsub = subset(GPSdf, DriftName == drift)
  NLsub =  subset(noiseDf, DriftName==drift)
  
  UTMflon <- approxfun(GPSsub$UTC, GPSsub$Longitude)
  UTMflat <- approxfun(GPSsub$UTC, GPSsub$Latitude)
  
  ##############################################################
  # Calculate the range from the whale to the GPS, TDOA and RL
  ############################################################
  
  # Lat/lon/ of the drift when the call was produced
  noiseDf$Lon[noiseDf$DriftName==drift] =UTMflon(NLsub$datetime_posix)
  noiseDf$Lat[noiseDf$DriftName==drift] = UTMflat(NLsub$datetime_posix)
  
  
}
colnames(noiseDf)[2]<-'NL'




ggplot(noiseDf)+
  geom_point(aes(x=Lon, y=Lat, color=NL))+
  scale_color_viridis_c(option='plasma')

# Clear out any NA values
noiseDf = noiseDf[!is.na(noiseDf$Lat),]

# add seconds since start as an integer for the gam
noiseDf$seconds = as.numeric(noiseDf$datetime_posix-median(noiseDf$datetime_posix))

library(fields)
loc<-cbind(noiseDf$Lon,noiseDf$Lat) #locations
obj<-spatialProcess(loc,noiseDf$NL, Distance = "rdist.earth",
                    cov.args = list(Covariance ="Exponential") )



###############################################################################
# Modelling Approach 1) GAMs
###############################################################################

# summarize by hour

HrlyNoise <- noiseDf %>%
  mutate(day_hour = format(datetime_posix, "%Y-%m-%d %H:00:00")) %>%  
  group_by(DriftName, day_hour) %>%  # Group by ID and day_hour
  summarise(NL = median(NL),
            Lat = median(Lat),
            Lon = median(Lon))  # Calculate median noise for each group

HrlyNoise$Date = as.POSIXct(result$day_hour, tz = 'UTC')

df <- df %>%
  left_join(median_hourly_noise, by = c("DriftName", "Hour"))


# add seconds since start as an integer for the gam
HrlyNoise$seconds = as.numeric(HrlyNoise$Date-median(HrlyNoise$Date))
# 
# # Step one, see how correlated the data are
# library(gratia)

mod1 = gam(data = result, NL ~ te(Lon, Lat)+ s(seconds))
mod2 = gam(data = result, NL ~ s(Lon)+s(Lat)+s(seconds))
mod3 = gam(data = result, NL ~ s(Lon)+s(Lat))

#draw(mod1)
#draw(mod2)


ilink_mod1 <- family(mod1)$linkinv
ilink_mod2 <- family(mod2)$linkinv
ilink_mod2 <- family(mod3)$linkinv


new_data=expand.grid(Lat = seq(min(GPSdf$Latitude)-.25, max(GPSdf$Latitude), 
                               length.out = 500),
                  Lon = seq(min(GPSdf$Longitude), max(GPSdf$Longitude)+.25, 
                            length.out = 500),
                  seconds = median(noiseDf$seconds))

#remove locations that were not within 10km of at least one drifter


new_data$MinDist = 0
for(ii in 1:nrow(new_data)){

  new_data$MinDist[ii] = min(haversine_dist(new_data$Lon[ii], new_data$Lat[ii],
                 GPSdf$Longitude, GPSdf$Latitude))
}

CombinedData =new_data[new_data$MinDist<= 10000,]
CombinedData$seconds=median(noiseDf$seconds)

pred <- predict(mod3, CombinedData, type = "link", se.fit = TRUE)
pred <- cbind(pred, CombinedData)
pred <- transform(pred, lwr_ci = ilink_mod3(fit - (2 * se.fit)),
                  upr_ci = ilink_mod1(fit + (2 * se.fit)),
                  fitted = ilink_mod1(fit))



library(ggplot2)
library(viridis)
ggplot(pred)+
  geom_tile(aes(x= Lon, y= Lat, fill = fit))+
  scale_fill_viridis_c(option='plasma', limits = c(82,90))+
  ggtitle('Ambient noie level mid-deployment')




CombinedData$seconds =quantile(noiseDf$seconds, .8)
CombinedData$seconds =median(noiseDf$seconds)


pred2 <- predict(mod2, CombinedData, type = "link", se.fit = TRUE)
pred2 <- cbind(pred2, CombinedData)
pred2 <- transform(pred, lwr_ci = ilink_mod2(fit - (2 * se.fit)),
                  upr_ci = ilink_mod2(fit + (2 * se.fit)),
                  fitted = ilink_mod2(fit))

ggplot(pred2)+
  geom_tile(aes(x= Lon, y= Lat, fill = fitted))+
  scale_fill_viridis_c(option='plasma', limits = c(82,90))+
  ggtitle('Ambient noie level Deployment End')


ggplot(result)+
  geom_point(aes(x=Lon, y=Lat, color=NL))+
  scale_color_viridis_c()



draw(mod1)

# Step two, make a tensor smooth plot of lat lon

