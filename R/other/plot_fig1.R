library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(tidyverse)


nc_data <- nc_open("C:/Users/a23092/Documents/Projects/Pelagic/IMRtag/data/current.nc")

print(nc_data)
lon <- ncvar_get(nc_data, "longitude")
length(lon)
lat <- ncvar_get(nc_data, "latitude")
length(lat)
time <- ncvar_get(nc_data, "time")
uo <- ncvar_get(nc_data, "uo")
dim(uo)

vo <- ncvar_get(nc_data, "vo")
dim(vo)

vo_df <- reshape2::melt(vo, varnames = c("lon", "lat", "time"), value.name = "NS_comp")
vo_df$lon <- as.vector(lon)[match(vo_df$lon, unique(vo_df$lon))]
vo_df$lat <- as.vector(lat)[match(vo_df$lat, unique(vo_df$lat))]
vo_df$time <- as.vector(time)[match(vo_df$time, unique(vo_df$time))]
uo_df <- reshape2::melt(uo, varnames = c("lon", "lat", "time"), value.name = "EW_comp")
uo_df$lon <- as.vector(lon)[match(uo_df$lon, unique(uo_df$lon))]
uo_df$lat <- as.vector(lat)[match(uo_df$lat, unique(uo_df$lat))]
uo_df$time <- as.vector(time)[match(uo_df$time, unique(uo_df$time))]

Dat <- vo_df
Dat$EW_comp <- uo_df$EW_comp
Dat$date <- as.POSIXct(Dat$time*3600, origin='1950-01-01 00:00:00', tz="GMT", format = "%Y-%m-%d")
Dat$speed <- sqrt(Dat$EW_comp^2 + Dat$NS_comp^2)

ggplot(Dat, aes(x=speed)) + geom_histogram()


Current_data1 <- Dat %>% filter(time %in% c(565008, 573768, 582552, 591312, 600072, 608832)) # June
Current_data1$Release_year <- (2014:2019)[match(Current_data1$time, c(565008, 573768, 582552, 591312, 600072, 608832))]
Current_data1 <- Current_data1 %>% filter(lat>51, lon<15, lat<68.5)
Current_data2 <- Dat %>% filter(time %in% c(564276, 573036, 581820, 590580, 599340, 608100)) # May
Current_data2$Release_year <- (2014:2019)[match(Current_data2$time, c(564276, 573036, 581820, 590580, 599340, 608100))]
Current_data2 <- Current_data2 %>% filter(lat>51, lon<15, lat<68.5)
Current_data3 <- Dat %>% filter(time %in% c(565740, 574500, 583284, 592044, 600804, 609564)) # July
Current_data3$Release_year <- (2014:2019)[match(Current_data3$time, c(565740, 574500, 583284, 592044, 600804, 609564))]
Current_data3 <- Current_data3 %>% filter(lat>51, lon<15, lat<68.5)
Current_data <- rbind(Current_data1, Current_data2)
Current_data <- Current_data %>% group_by(Release_year,lat,lon) %>% summarise(EW_comp = mean(EW_comp), NS_comp = mean(NS_comp))
Current_data$speed <- sqrt(Current_data$EW_comp^2 + Current_data$NS_comp^2)

#Current_data <- Current_data1

library(sf)
world <- st_read("C:/Users/a23092/Documents/Projects/Mackerel_distribution/Shapefiles/ne_10m_land.shp")
Atlantic <- st_crop(world, c(xmin = -51.5, ymin = 50, xmax = 32, ymax = 80))

#2014
g2014_05 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==564276, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
g2014_06 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==565008, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
ggsave(g2014_05, file=paste0(getwd(), "/plots/current_map/2014_05.pdf"), width=14, height = 14, dpi=450)
ggsave(g2014_06, file=paste0(getwd(), "/plots/current_map/2014_06.pdf"), width=14, height = 14, dpi=450)
#2015
g2015_05 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==573036, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
g2015_06 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==573768, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
ggsave(g2015_05, file=paste0(getwd(), "/plots/current_map/2015_05.pdf"), width=14, height = 14, dpi=450)
#2016
g2016_05 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==581820, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
g2016_06 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==582552, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
ggsave(g2016_05, file=paste0(getwd(), "/plots/current_map/2016_05.pdf"), width=14, height = 14, dpi=450)
#2017
g2017_05 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==590580, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
g2017_06 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==591312, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
ggsave(g2017_05, file=paste0(getwd(), "/plots/current_map/2017_05.pdf"), width=14, height = 14, dpi=450)
#2018
g2018_05 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==599340, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
g2018_06 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==600072, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
ggsave(g2018_05, file=paste0(getwd(), "/plots/current_map/2018_05.pdf"), width=14, height = 14, dpi=450)
#2019
g2019_05 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==608100, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
g2019_06 <- ggplot(Atlantic) + geom_sf() + theme_bw() +
  geom_segment(data = Dat %>% filter(time==608832, speed>0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.05, "inches"), type = "closed"), alpha = 0.4) + scale_color_gradient(low="lightblue1", high="darkblue")
ggsave(g2019_05, file=paste0(getwd(), "/plots/current_map/2019_05.pdf"), width=14, height = 14, dpi=450)

