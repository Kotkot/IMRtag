# ---------------------------------
# ----------------------------------
# KOTARO:
#  To run this script, remember to first run main.R until about line 170
#  The Figure 1 is created from here to ~ line 115
# ----------------------------------

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(tidyverse)
library(ggnewscale)


nc_data <- nc_open("data/current.nc")
nc_data <- nc_open("data/current_new.nc")

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

### May current data
Current_data <- Dat %>% filter(time %in% c(564276, 573036, 581820, 590580, 599340, 608100, 616884))
Current_data$Release_year <- (2014:2020)[match(Current_data$time, c(564276, 573036, 581820, 590580, 599340, 608100, 616884))]
Current_data <- Current_data %>% filter(lat>51, lon<15, lat<68.5)


lag = 0
years <- 2014:2020
library(lubridate)
library(sf)
dat_release <- Data_mackerel_all %>% subset(Tag_area %in% c("South_Ireland", "North_Ireland"))
dat_release$julian <-  sapply(1:nrow(dat_release), function(k) as.numeric(julian(dat_release$ReleaseDate[k],
                                                                                 as.POSIXct(paste0(dat_release$Release_year[k], "-01-01"), tz = "GMT"))))
dat_release$Tag_area <- as.character(dat_release$Tag_area)
dat_release$Tag_area <- as.factor(dat_release$Tag_area)
dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))

dat_recap <- Data_mackerel_final
dat_recap$Release_year <- as.integer(as.character(dat_recap$Release_year))
dat_recap <- dat_recap %>% subset(Catch_year==Release_year+lag &
                                              Tag_area %in% c("South_Ireland", "North_Ireland"))
dat_recap$julian <-  sapply(1:nrow(dat_recap), function(k) as.numeric(julian(dat_recap$RecaptureDate[k], as.POSIXct(paste0(dat_recap$Release_year[k], "-01-01"), tz = "GMT"))))
dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)

dat_release <- filter(dat_release, Release_year %in% years)
dat_recap <- filter(dat_recap, Release_year %in% years)

dat_recap_new <- as_tibble(dat_recap)%>% dplyr::select(cLon, cLat, julian, Release_year) %>%
  group_by(Release_year) %>% count(round(cLon,1), round(cLat,1), round(julian,0))
colnames(dat_recap_new) <- c("Release_year", "cLon", "cLat", "julian_recapture", "n")




dat_release_new <- as_tibble(dat_release) %>% dplyr::select(Longitude, Latitude, julian, Release_year) %>%
  group_by(Release_year) %>% count(round(Longitude,1), round(Latitude,1), round(julian,0))
colnames(dat_release_new) <- c("Release_year", "Longitude", "Latitude", "julian_release", "n")


ncount_release <- dat_release %>% group_by(Release_year) %>% count()
ncount_release$lon = c(-32)
ncount_release$lat = c(53)
ncount_release$label = paste0("N-released = ", ncount_release$n)
ncount_recap <- dat_recap %>% group_by(Release_year) %>% count()
ncount_recap$lon = c(-16)
ncount_recap$lat = c(70)
ncount_recap$label = paste0("N-recaptured = ", ncount_recap$n, "\n (rate = ", round(ncount_recap$n/ncount_release$n,3)*100, "%)")
Select_catch <- filter(Catch_data, year(CatchDate) %in% as.integer(years+lag))
Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
Select_catch <- subset(Select_catch, cLat<72)
# --- Should select take out catches prior <= may   ---
# --- but this info is not in Catch_data            ---
#Select_catch <- filter(Select_catch) # catches from june - december
hull_catch <- NULL
for(i in 1:length(years)){
(hull_catch <- rbind(hull_catch, filter(Select_catch, year(CatchDate) == years[i],
                                        month(ProcessingDate) %in% 6:12) %>%
                       dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat, CatchDate) %>%
  mutate(Release_year = year(CatchDate))))
}


(weights <- Catch_data %>% filter(between(month(CatchDate), 6,12)) %>% mutate(Catch_year = year(CatchDate)) %>%
  group_by(Catch_year) %>% summarize(biomass = sum(CatchWeight/1e6)) %>% filter(between(Catch_year,2014,2020)) %>%
  rename(Release_year = Catch_year) %>% add_column(xpos = -Inf, ypos = Inf)
)
area <- Polygon(hull_catch[,c('cLon','cLat')])@area
Norway <- st_read("shapefile/ne_10m_land.shp")
Norway <- st_crop(Norway, c(xmin = -31, ymin = 50, xmax = 19, ymax = 70))

#dat_release_new[which.max(dat_release_new$lo) ,"lo"] <- -dat_release_new[which.max(dat_release_new$lo) ,"lo"]
g1 <- ggplot(Norway) + geom_sf() +
  # geom_segment(data = Current_data %>% filter(speed > 0.1), aes(x=lon, y=lat, xend = lon + EW_comp*5, yend = lat + NS_comp*5, col=speed),
  #              arrow = arrow(angle = 15, length = unit(0.03, "inches"), type = "closed"), alpha = 0.8, size=0.3) +
  # scale_color_gradient(low="lightblue1", high="darkblue",  name = "Current\nspeed (m/s)") +
  # new_scale_color() +
  # geom_point(data=dat_release_new %>% arrange(desc(n)),
  #            aes(x=Longitude, y=Latitude,col=julian_release, size=n)) +
  geom_point(data=dat_release_new %>% arrange(desc(n)),
             aes(x=Longitude, y=Latitude,col=julian_release),size=2) +
  #geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat)) +
  scale_color_viridis_c(limits = c(120, 160), name = "Julian\nRelease",
                        breaks = seq(120,160,5),
                        guide = guide_colorbar(barheight = 15, ticks.colour = "grey20")) +
  #guides(color = guide_colorbar(barheight = 15))+
  # scale_size_continuous(breaks=c(2,5,10,100,1000,5000),range = c(1, 7), name = "#Releases", trans="log") +
  scale_size_continuous(breaks=seq(1,7),range = c(1, 7), name = "#Recaptures") +
  new_scale_color() +
  geom_point(data=dat_recap_new %>% arrange(desc(julian_recapture)),
             aes(x=cLon, y=cLat, col = julian_recapture, size=n),
                 #size = 900+200*order(julian_recapture, decreasing = FALSE)/nrow(dat_recap_new)),
             pch = 18) +
  scale_color_gradient2(low ="darkgreen", mid = "yellow", high = "red", name = "Julian\nRecapture",
                        limits = c(180,300),
                        midpoint = (300+180)/2,
                        na.value = "red",
                        guide = guide_colorbar(barheight = 15, ticks.colour = "grey20"),
                        breaks = seq(180,300,20),
                        minor_breaks = seq(180,300,10)) +
  #scale_size_continuous(range = c(0.3, 0.6)) +
  scale_x_continuous(expand =c(0,0))+scale_y_continuous(expand = c(0,0))+
  guides(size = guide_legend(override.aes = list(shape = 19)))+
  #geom_text(data=data.frame(x=-Inf,y=60,label="N-recapture = "), aes(x=x, y=y, label=label), col="black",
  #          fontface="bold") +
  #geom_text(data=data.frame(x=-26,y=54,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black",
  #          fontface="bold") +
  geom_text(data =weights, aes(x = xpos, y = ypos, label = paste0("Biomass scanned: ", round(biomass,1), " kt")),
            vjust = 1.2, hjust = -.75, fontface = "bold", col = "black")+
  geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA, lty  =2) +
  facet_wrap(~Release_year, ncol = 2)+
  geom_text(data=ncount_release, aes(x=5, y=52.5, label=label), hjust=-.05, size=3)+
  geom_text(data=ncount_recap, aes(x=5, y=51, label=label), hjust=-.05, size=3) + labs(x="Longitude", y="Latitude") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_text(face = "bold", size = 14),
                     axis.text = element_text(size=12),
                     axis.title = element_text(size=16),
                     strip.background = element_rect(fill = "white"),
                     strip.text = element_text(face ="bold", size = 16)) +
  geom_hline(yintercept=62, linetype=2) + geom_vline(xintercept=-10, linetype=2)
                     #panel.spacing.x = unit(6.5, "mm"))+
ggsave(g1, filename = "MS/figs/Fig1_May_average_new.pdf",
       width=42, height=32, units="cm", dpi = 450)
ggsave(g1, filename = "MS/figs/Fig1_May_average_new1.pdf",
       width=32, height=42, units="cm", dpi = 450)


##### Figure 2 in the paper
Catch_data_test <- Catch_data %>% mutate(date = as.Date(ProcessingDate))
Catch_data_test$Year <- as.numeric(substr(Catch_data_test$date, start=1, stop=4))
Catch_data_test <- Catch_data_test %>% filter(Year>=2014)
Catch_data_test$julian_catch <-  as.numeric(julian(as.POSIXct(Catch_data_test$date), as.POSIXct(paste0(2014, "-01-01"), tz = "GMT")))
Catch_data_test$julian_catch_std <-  Catch_data_test$julian_catch %% 365
Catch_data_test$month <-  month(as.POSIXct(Catch_data_test$date))
Catch_data_test <- Catch_data_test %>% filter(julian_catch_std >= 180) %>% group_by(Year, ICES_Rectangle) %>% summarise(Catch=sum(CatchWeight), Julian = mean(julian_catch_std))
library(mapplots)
coords = ices.rect(Catch_data_test$ICES_Rectangle)
Catch_data_test <- cbind(Catch_data_test, coords)
Catch_data_test <- Catch_data_test %>% mutate(Catch_ton = Catch/1000)


plot_catch <- ggplot(Norway) + geom_sf() +
  geom_tile(data=Catch_data_test, aes(x=lon, y=lat, fill=Catch_ton)) +
  facet_wrap(~Year, ncol=2) + scale_fill_gradient(low="lightblue1", high="darkblue", name = "Catch (t)") +
  labs(x="Longitude", y="Latitude") +
  coord_sf(xlim=c(-30,20), ylim=c(50,70)) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_text(face = "bold", size = 13),
                     strip.background = element_rect(fill = "white"),
                     strip.text = element_text(face ="bold", size = 15),
                     axis.title = element_text(size=15),
                     axis.text = element_text(size=13))

plot_movement <- ggplot(Norway) + geom_sf() +
  geom_point(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat, col=log_rate)) +
  facet_wrap(~Release_year, ncol=2) + scale_color_gradient(low="lightblue1", high="darkblue", name = "Mvt rate\n   (log)") +
  labs(x="Longitude", y="Latitude") +
  coord_sf(xlim=c(-30,20), ylim=c(50,70)) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_text(face = "bold", size = 13),
                     strip.background = element_rect(fill = "white"),
                     strip.text = element_text(face ="bold", size = 15),
                     axis.title = element_text(size=15),
                     axis.text = element_text(size=13))

Fig2 <- grid.arrange(plot_catch, plot_movement, ncol=2)

ggsave(Fig2, filename = "MS/figs/Fig2.pdf",
       width=15, height=12, dpi = 450)







# -----------------------------------------------------------------------
# ---- Another figure with the location of the Channel without the tags
# -----------------------------------------------------------------------

library(ggOceanMaps)

Norway <- st_read("shapefile/ne_10m_land.shp")
Norway <- st_crop(Norway, c(xmin = -34, ymin = 49, xmax = 22, ymax = 72))

dt <- data.frame(lon = c(-31, -31, 19, 19), lat = c(45, 70, 70, 45))

Current_data1 <- Dat %>% filter(time %in% c(565008, 573768, 582552, 591312, 600072, 608832)) # June
Current_data1$Release_year <- (2014:2019)[match(Current_data1$time, c(565008, 573768, 582552, 591312, 600072, 608832))]
Current_data1 <- Current_data1 %>% filter(lat>51, lon<15, lat<68.5)
Current_data2 <- Dat %>% filter(time %in% c(564276, 573036, 581820, 590580, 599340, 608100)) # May
Current_data2$Release_year <- (2014:2019)[match(Current_data2$time, c(564276, 573036, 581820, 590580, 599340, 608100))]
Current_data2 <- Current_data2 %>% filter(lat>51, lon<15, lat<68.5)
Current_data3 <- Dat %>% filter(time %in% c(565740, 574500, 583284, 592044, 600804, 609564)) # July
Current_data3$Release_year <- (2014:2019)[match(Current_data3$time, c(565740, 574500, 583284, 592044, 600804, 609564))]
Current_data3 <- Current_data3 %>% filter(lat>51, lon<15, lat<68.5)
Current_data <- Current_data2

survey <- Current_data %>% st_as_sf(crs = 4326, coords = c("lon", "lat"))
surv_utm_coords <- st_coordinates(survey)
Current_data$X <- surv_utm_coords[,1]
Current_data$Y <- surv_utm_coords[,2]

X <- basemap_data(limits = NULL, data = dt, shapefiles = NULL,
                  bathymetry = TRUE, glaciers = FALSE, resolution = "low",
                  lon.interval = NULL, lat.interval = NULL,
                  rotate = FALSE, verbose = FALSE)

Bathy <- X$shapefiles$bathy %>% st_as_sf(crs = 4326, coords = c("lon", "lat"))

cc <- scales::seq_gradient_pal("lightblue1", "darkblue", "Lab")(seq(0,1,length.out=8))

g2 <- ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
  geom_sf(data = Bathy, aes(col=depth, fill=depth)) +
  scale_colour_manual(values=cc, name = "Depth") + scale_fill_manual(values=cc, name = "Depth") + geom_sf(data=Norway, col=grey(0.3)) +
  coord_sf(xlim=c(-31, 19), ylim=c(50, 70)) +
  new_scale_color() +
  new_scale_fill() +
  geom_segment(data = Current_data %>% filter(speed > 0.1), aes(x=X, y=Y, xend = X + EW_comp*5, yend = Y + NS_comp*5, col=speed),
               arrow = arrow(angle = 15, length = unit(0.03, "inches"), type = "closed"), alpha = 0.8, size=0.3) +
  scale_color_gradient(low="lightpink1", high="darkred", name = "Current\nspeed (m/s)") +
  facet_wrap(~Release_year, ncol = 2) + labs(x="Longitude", y="Latitude") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_text(face = "bold", size = 12),
                     strip.background = element_rect(fill = "white"),
                     strip.text = element_text(face ="bold", size = 14))
#panel.spacing.x = unit(6.5, "mm"))+
ggsave(g2, filename = "MS/figs/Fig2_May_Current_bathym.pdf",
       width=28, height=32, units="cm", dpi = 450)



g3 <- ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
  geom_sf(data = Bathy, aes(col=depth, fill=depth)) + theme_bw() +
  scale_colour_manual(values=cc) + scale_fill_manual(values=cc) + geom_sf(data=Norway, col=grey(0.3)) +
  coord_sf(xlim=c(-31, 19), ylim=c(50, 70)) + labs(x="Longitude", y="Latitude") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_text(face = "bold", size = 12),
                     strip.background = element_rect(fill = "white"),
                     strip.text = element_text(face ="bold", size = 14))
#panel.spacing.x = unit(6.5, "mm"))+
ggsave(g3, filename = "MS/figs/Fig_background.pdf",
       width=18, height=18, units="cm", dpi = 450)











# ------------------------------------------
# ---- Arils supplementary figure ----------
# ------------------------------------------
dat_release <- Data_mackerel_all %>% subset(Tag_area %in% c("South_Ireland", "North_Ireland"))
dat_release$julian <-  sapply(1:nrow(dat_release), function(k) as.numeric(julian(dat_release$relesedate[k],
                                                                                 as.POSIXct(paste0(dat_release$Release_year[k], "-01-01"), tz = "GMT"))))
dat_release$Tag_area <- as.character(dat_release$Tag_area)
dat_release$Tag_area <- as.factor(dat_release$Tag_area)
dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))


ncount_release <- filter(dat_release, Release_year %in% years) %>% group_by(Release_year, Tag_area, Release_timing) %>% count() %>%
  mutate(latitudes = ifelse(Tag_area == "From N. Ireland", ">54°N", "<54°N"),
         timing = ifelse(Release_timing == "First_half", "<141", ">= 141"))
ncount_release$lon = -Inf
ncount_release$lat = ifelse(ncount_release$Tag_area == "From N. Ireland" & ncount_release$timing == "<141", 62,
                            ifelse(ncount_release$Tag_area != "From N. Ireland" & ncount_release$timing == "<141",61,
                                   ifelse(ncount_release$Tag_area == "From N. Ireland" & ncount_release$timing != "<141",60,
                                          59)))
ncount_release <- add_column(ncount_release,
                             label =  paste0("N-released lat", ncount_release$latitudes, ", julian day",
                                             ncount_release$timing, "= ", ncount_release$n))


ncount_release <- filter(dat_release, Release_year %in% years) %>%
  group_by(Release_year, Release_timing) %>% count() %>%
  mutate(timing = ifelse(Release_timing == "First_half", "<141", ">=141"))
ncount_release_julian$Release_timing = NULL
ncount_release_julian <- spread(ncount_release_julian, timing, n)
ncount_release_julian$label  <-NA
for(i in 1:nrow(ncount_release_julian)){
  ncount_release_julian$label[i] <- paste0("N-released julian day  <141 = ", ncount_release_julian[i,'<141'],
                                  "\nN-released julian day >=141 = ", ncount_release_julian[i,'>=141'])
}
ncount_release_julian$lon <- -Inf
ncount_release_julian$lat <- 60

ncount_release_NS <- filter(dat_release, Release_year %in% years) %>%
  group_by(Release_year, Tag_area) %>% count() %>%
  mutate(latitudes = ifelse(Tag_area == "From N. Ireland", ">54°N", "<54°N"))
ncount_release_NS$Tag_area = NULL
ncount_release_NS <- spread(ncount_release_NS, latitudes, n)
ncount_release_NS$label  <-NA
for(i in 1:nrow(ncount_release_NS))
ncount_release_NS$label  <- paste0("N-released lat >54°N = ", ncount_release_NS[i,'>54°N'],
                                   "\nN-released lat>=54°N = ", ncount_release_NS[i,'<54°N'])
ncount_release_NS$lon <- -Inf
ncount_release_NS$lat <- 58




dat_recap <- Data_mackerel_final
dat_recap$Release_year <- as.integer(as.character(dat_recap$Release_year))
dat_recap <- dat_recap %>% subset(Catch_year==Release_year+lag &
                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
dat_recap$julian_release_per_year <-  sapply(1:nrow(dat_recap), function(k) as.numeric(julian(dat_recap$relesedate[k], as.POSIXct(paste0(dat_recap$Release_year[k], "-01-01"), tz = "GMT"))))
dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)

dat_release <- filter(dat_release, Release_year %in% years)
dat_recap <- filter(dat_recap, Release_year %in% years)

dat_recap_new <- as_tibble(dat_recap)%>% dplyr::select(cLon, cLat, julian_release_per_year, Release_year,Tag_area) #%>%
  #group_by(Release_year,Tag_area) # %>%count(round(cLon,1), round(cLat,1), round(julian_release_per_year,0)) %>%
dat_recap_new$Tag_area <-   mapvalues(dat_recap_new$Tag_area, from  = c("South_Ireland", "North_Ireland"),
                  to = c("South Ireland", "North Ireland"))
#colnames(dat_recap_new) <- c("Release_year", "Tag_area",  "julian_release")#, "n")

Norway <- st_read("shapefile/ne_10m_land.shp")
Norway <- st_crop(Norway, c(xmin = min(dat_recap_new$cLon)-1,
                            ymin = min(dat_recap_new$cLat)-1,
                            xmax = max(dat_recap_new$cLon)+3,
                            ymax = max(dat_recap_new$cLat)+1))
range(dat_recap_new$cLon)
set.seed(13112020)

dat_recap_new <- left_join(dat_recap_new, dat_recap_new %>% group_by(cLon,cLat, Release_year) %>% count(), by = c("cLon","cLat", "Release_year"))
dat_recap_new$width <- dat_recap_new$height <- ifelse(dat_recap_new$n == 1, 0,
                                                      ifelse(dat_recap_new$n < 10, .3,
                                                             ifelse(dat_recap_new$n <20, .4, .5)))

dat_recap_new$cLonjit <- jitter(dat_recap_new$cLon, amount = resolution(dat_recap_new$cLon, zero = FALSE) * dat_recap_new$width)
dat_recap_new$cLatjit <- jitter(dat_recap_new$cLat, amount = resolution(dat_recap_new$cLat, zero = FALSE) * dat_recap_new$height)
dat_recap_new$Tag_area <- factor(dat_recap_new$Tag_area, labels=c(">=54°N", "<54°N"))
g2 <- ggplot(Norway) + geom_sf() +
  geom_point(data=dat_recap_new %>% arrange((julian_release_per_year)), aes(x = cLonjit, y = cLatjit, col = julian_release_per_year,
                                     shape = Tag_area),
             size = 2.5)+
             #position = "jitter",
             #width = 0.3,
             #height =0.3) +

  #geom_text(data = ncount_release_julian, aes(x = lon, y = lat-.4, label = label),
  #          hjust = 0)+
 # geom_text(data = ncount_release_NS, aes(x = lon, y = lat, label = label),
  #          hjust = 0)+
  geom_text(data = data.frame(x = c(-20,2,2),
                              y = c(67.5, 67.5, 58),
                              lab = c("Norwegian Sea\nWest",
                                      "Norwegian Sea\nEast",
                                      "North Sea"),
                              Release_year = rep(2014,3)),
            aes(x=x, y =y, label = lab), color = "blue",
            size = 5) +
  facet_wrap(~Release_year, ncol = 2)+
    scale_color_viridis_c(name = "Julian day\nof release", breaks = c(seq(120,160,5),141)) +
    scale_x_continuous(expand = c(0,0), name = "Longitude") +
    scale_y_continuous(expand = c(0,0), name = "Latitude") +
  scale_shape_discrete(name = "Tagging area", solid = TRUE) +
    guides(color = guide_colorbar(barheight = 35)) +
  geom_hline(aes(yintercept = 62), lty = 2, col = "blue")+
  geom_vline(aes(xintercept = -5), lty = 2, col = "blue")+
  theme_bw()+
  theme(strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(face = "bold", size = 12))

#g2
ggsave(g2, filename = "plots/Arils_supplementary_figure_jitter_without_text.tiff",
       width=28, height=30, units="cm", device = "tiff", dpi = "retina")


# ---------------------------------------------------------------------------------
# Same figure but make recaptures based on releases split into 5 body length bins,
# 5 different size of symbols, biggest for biggest length bin, and make sure biggest
# always are on bottom smallest on top at each position (jittering not necessary) to
# demonstrate if the length bin was recaptured in that position or not


quantile(dat_recap$length,seq(.2,.8,.2))
table(cut_number(dat_recap$length, n = 5))
table(cut(dat_recap$length, breaks = c(-Inf,33,35,37,39, Inf)))
# Update length bin figure with 33 and smaller, 34-35, 36-37, 38-39, 40 and bigger
dat_recap$length_cat <- cut(dat_recap$length, breaks = c(min(dat_recap$length),33,35,37,39, max(dat_recap$length)),
                            include.lowest = T)
dat_recap$length_int <- as.integer(dat_recap$length_cat)
dat_recap$Tag_area <-   plyr::mapvalues(dat_recap$Tag_area, from  = c("South_Ireland", "North_Ireland"),
                                      to = c("<54°N", ">54°N"))

g3 <- ggplot(Norway) + geom_sf() +
  geom_point(data=dat_recap%>% arrange(julian_release_per_year, length_cat),
             aes(x = cLon, y = cLat, fill = julian_release_per_year,
                 #shape = Tag_area,
                 size = length_cat),
             pch = 21)+
  #position = "jitter",
  #width = 0.3,
  #height =0.3) +

  # geom_text(data = ncount_release_julian, aes(x = lon, y = lat-.4, label = label),
  #           hjust = 0)+
  # geom_text(data = ncount_release_NS, aes(x = lon, y = lat, label = label),
  #           hjust = 0)+
  facet_wrap(~Release_year, ncol = 2)+
  scale_fill_viridis_c(name = "Julian day\nof release", breaks = c(seq(120,160,5))) +
  scale_x_continuous(expand = c(0,0), name = "Longitude") +
  scale_y_continuous(expand = c(0,0), name = "Latitude") +
  scale_size_discrete(name = "Length bin")+
  #scale_shape_discrete(solid = FALSE) +
  scale_shape_discrete(name = "Released", solid = TRUE) +
  guides(fill = guide_colorbar(barheight = 35)) +
  theme_bw()+
  theme(strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(face = "bold", size = 12))
ggsave(g3, filename = "plots/Arils_supplementary_figure_length_bins.tiff",
       width=28, height=30, units="cm", device = "tiff", dpi = "retina")

dat_release$length_cat <- cut(dat_release$length, breaks = c(min(dat_release$length),33,35,37,39, max(dat_release$length)),
                              include.lowest =T)
ggplot(dat_release %>%
         mutate(Tag_area = ifelse(Tag_area == "From N. Ireland", ">54°N", "<54°N")), aes(x=factor(length))) + geom_bar(aes(fill = Tag_area)) +
  facet_wrap(~Release_year, ncol = 1, strip.position = "right")+
  theme_bw() +
  scale_fill_discrete(name = "Tagging area") +
  xlab ("Length") + ylab("Frequency")+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "top",
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(face = "bold", size = 12))#+ scale_y_continuous(expand = c(0,0))
  #stat_density(alpha = 0.1, fill = "red", col = "red")
ggsave("plots/length_cat_distribution_releases.tiff",
       width=25, height=30, units="cm", device = "tiff", dpi = "retina")

# -- multipanel figure
# Upper left: numbers released per year stacked with different colors for numbers north-south, upper right same but stacked with different colours for the before and after 22May
# Mid left: numbers recaptured per year and area (3 columns per year) stacked with different colors for numbers north-south, mid right same but stacked with different colours for the before and after 22May
# Bottom left: rate (numbers recaptured/released) per year and area stacked with different colors for rate releasaed north-south, bottom right same but stacked with different colours for the before and after 22May
# -----------------

dat_release <- dat_release %>% mutate(recaptureArea =
                                        ifelse(cLat <= 62, "NS",
                                               ifelse(cLon>5, "E", "W")))

ncount_release <- filter(dat_release, Release_year %in% years) %>%
  group_by(Release_year, Tag_area, Release_timing,recaptureArea) %>% count() %>%
  mutate(latitudes = ifelse(Tag_area == "From N. Ireland", ">=54°N", "<54°N"),
         timing = ifelse(Release_timing == "First_half", "<141", ">= 141"))
ncount_release$lon = -Inf
ncount_release$lat = ifelse(ncount_release$Tag_area == "From N. Ireland" & ncount_release$timing == "<141", 62,
                            ifelse(ncount_release$Tag_area != "From N. Ireland" & ncount_release$timing == "<141",61,
                                   ifelse(ncount_release$Tag_area == "From N. Ireland" & ncount_release$timing != "<141",60,
                                          59)))
ncount_release <- add_column(ncount_release,
                             label =  paste0("N-released lat", ncount_release$latitudes, ", julian day",
                                             ncount_release$timing, "= ", ncount_release$n))
# -------------------------------------------------
# --------------- Figure top left :  --------------
# - Numbers released per before and after May 22 --
# -------------------------------------------------
ncount_release_julian <- filter(dat_release, Release_year %in% years) %>%
  dplyr::select(Release_year, Release_timing) %>%
  group_by(Release_year, Release_timing) %>% count() %>%
  mutate(timing = ifelse(Release_timing == "First_half", "<141", ">=141")) %>%
  dplyr::select(-Release_timing) %>% rename(freq = n)

ncount_release_julian$totFreq <- NA
for(i in unique(ncount_release_julian$Release_year)){
  ncount_release_julian$totFreq[ncount_release_julian$Release_year == i &
                                  ncount_release_julian$timing == '>=141'] <- ncount_release_julian$freq[ncount_release_julian$Release_year == i &

                                                                                                          ncount_release_julian$timing == '<141'] +
    ncount_release_julian$freq[ncount_release_julian$Release_year == i &

                                 ncount_release_julian$timing == '>=141']/2
  ncount_release_julian$totFreq[ncount_release_julian$Release_year == i &
                                  ncount_release_julian$timing == '<141'] <- ncount_release_julian$freq[ncount_release_julian$Release_year == i &
                                                                                                              ncount_release_julian$timing == '<141'] /2
}

# turning the order:
ncount_release_julian$timing <- factor(ncount_release_julian$timing, levels = c(">=141" = ">=141", "<141"="<141"))

(p1 <- ggplot(ncount_release_julian, aes(x = Release_year, fill = timing, y = freq)) +
    geom_bar(position="stack", stat="identity", width = .99) +
    geom_text(aes(label=freq, y = totFreq))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top",
          axis.text.y = element_text(angle = 90, hjust = .5),
          axis.text = element_text(size = 11),
          #strip.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size =14),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())+
    scale_x_discrete(expand = c(0,0), name = "Release year")+
    scale_y_continuous(expand = c(0,500), "Number of released fish") +
    scale_fill_discrete("Julian day of release"))

# -------------------------------------------------
# --------------- Figure top right :  -------------
# - Numbers released per north or south of 54N ----
# -------------------------------------------------

ncount_release_NS <- filter(dat_release, Release_year %in% years) %>%
  dplyr::select(Release_year, Tag_area) %>%
  group_by(Release_year, Tag_area) %>% count() %>%
  mutate(timing = ifelse(Tag_area == "From N. Ireland", ">=54°N", "<54°N")) %>%
  dplyr::select(-Tag_area) %>% rename(freq = n) %>% mutate(freq = freq)

ncount_release_NS$totFreq <- NA
for(i in unique(ncount_release_NS$Release_year)){
  ncount_release_NS$totFreq[ncount_release_NS$Release_year == i &
                                  ncount_release_NS$timing != "<54°N"] <- ncount_release_NS$freq[ncount_release_NS$Release_year == i &

                                                                                                          ncount_release_NS$timing == "<54°N"] +
    ncount_release_NS$freq[ncount_release_NS$Release_year == i &

                                 ncount_release_NS$timing != "<54°N"]/2
  ncount_release_NS$totFreq[ncount_release_NS$Release_year == i &
                                  ncount_release_NS$timing == "<54°N"] <- ncount_release_NS$freq[ncount_release_NS$Release_year == i &
                                                                                                           ncount_release_NS$timing == "<54°N"] /2
}
ncount_release_NS <- ncount_release_NS %>% mutate(Tag_area = timing, timing = NULL)
ncount_release_NS$Tag_area <- factor(ncount_release_NS$Tag_area,
                                     levels = c(">=54°N" = ">=54°N", "<54°N"="<54°N"))

(p2 <- ggplot(ncount_release_NS, aes(x = Release_year, fill = Tag_area, y = freq)) +
    geom_bar(position="stack", stat="identity", width = .99) +
    geom_text(aes(label=freq, y = totFreq))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top",
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(size =14),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())+
    scale_x_discrete(expand = c(0,0), name = "Release year")+
    scale_y_continuous(expand = c(0,500), "Number of released fish") +
    scale_fill_discrete("Tagging area"))

# --------------------------------------------------
# --------------- Figure mid left :  ---------------
# - Numbers recaptured per before and after May 22 -
# --------------------------------------------------
dat_recap<- dat_recap %>% mutate(recaptureArea =
                                        ifelse(cLon < -5 & cLat >= 62, "W",
                                               ifelse(cLon>=-5 & cLat >= 62, "E",
                                                      ifelse(cLon >=-5 & cLat < 62, "NS", "Other"))))
ncount_recap_julian <- filter(dat_recap, Release_year %in% years) %>%
  dplyr::select(Release_year, Release_timing, recaptureArea) %>%
  group_by(Release_year, Release_timing,recaptureArea) %>% count() %>%
  mutate(timing = ifelse(Release_timing == "First_half", "<141", ">=141")) %>%
  dplyr::select(-Release_timing) %>% rename(freq = n)
ncount_recap_julian$Release_timing <- NULL
ncount_recap_julian <-ncount_recap_julian %>% group_by(Release_year, recaptureArea) %>%
  mutate(relFreq = freq/sum(freq), totFreqRecap = sum(freq))


  ncount_recap_julian$yeararea = factor(str_replace(as.character(interaction(ncount_recap_julian$Release_year,
                                                              ncount_recap_julian$recaptureArea)),
                                     '\\.', ' / '), ordered=TRUE)

(p3 <- ggplot(ncount_recap_julian, aes(x =yeararea, fill = timing, y = freq)) +
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label=freq, y = totFreq))+
    #facet_wrap( ~ Release_year, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(face = "bold"),
          #panel.margin = grid::unit(-1.25, "lines"),
          #strip.placement = "outside",
          #strip.background = element_rect(fill = "white", color = "white"),
          axis.text.x = element_text(angle = 90),
          strip.text = element_blank())+
    scale_x_discrete(expand = c(0,0), name = "Release year/area of recapture")+
    scale_y_continuous(expand = c(0,2), "Number of recaptured fish") +
    scale_fill_discrete("Julian\nRelease"))

  ncount_recap_julian <- transform(ncount_recap_julian, recaptureArea = factor(recaptureArea, levels = c("W","E","NS")),
                                   timing = factor(timing, levels = c(">=141","<141")))

(p3 <- ggplot(ncount_recap_julian, aes(x = recaptureArea, fill = timing, y = freq)) +
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label=freq, y = totFreq))+
    facet_wrap( ~ Release_year, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(face = "bold"),
          panel.margin = grid::unit(0, "lines"),
          strip.placement = "outside",
          legend.position = "none",
          #axis.text = element_text(size = 12),
          #strip.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(angle = 90, hjust = .5, size = 11),
          strip.background = element_rect(fill = "white", color = "white"),
          #axis.text.x = element_text(angle = 0),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.text = element_blank()
          )+
    scale_x_discrete(expand = c(.2,.2), name = "Release year")+
    scale_y_continuous(expand = c(0,0), "Number of recaptured fish", breaks = seq(0,200,25)) +
    scale_fill_discrete("Julian\n release"))
# --------------------

# --------------------------------------------------
# --------------- Figure mid left :  ---------------
# - Numbers recaptured per before and after May 22 -
# --------------------------------------------------

  levels(dat_recap$Tag_area ) <- c(">=54°N","<54°N")
ncount_recap_NS <- filter(dat_recap, Release_year %in% years) %>%
  dplyr::select(Release_year, Tag_area, recaptureArea) %>%
  group_by(Release_year, Tag_area,recaptureArea) %>% count()  %>% rename(freq = n)#%>%

  ncount_recap_NS <- ncount_recap_NS %>%  group_by(Release_year, recaptureArea) %>%
    mutate(relFreq = freq/sum(freq), totFreqRecap = sum(freq))



  #mutate(timing = ifelse(Tag_area ==  "From N. Ireland", ">54°N", "<54°N")) %>%
  #dplyr::select(-Tag_area)

# ncount_recap_NS$totFreq <- NA
# for(i in unique(ncount_recap_NS$Release_year)){
#   ncount_recap_NS$totFreq[ncount_recap_NS$Release_year == i &
#                                 ncount_recap_NS$timing == '<141'] <- ncount_recap_NS$freq[ncount_recap_NS$Release_year == i &
#
#                                                                                                     ncount_recap_NS$timing == '>=141'] +
#     ncount_recap_NS$freq[ncount_recap_NS$Release_year == i &
#
#                                ncount_recap_NS$timing == '<141']/2
#   ncount_recap_NS$totFreq[ncount_recap_NS$Release_year == i &
#                                 ncount_recap_NS$timing == '>=141'] <- ncount_recap_NS$freq[ncount_recap_NS$Release_year == i &
#                                                                                                      ncount_recap_NS$timing == '>=141'] /2
# }

ncount_recap_NS$yeararea <- factor(str_replace(interaction(ncount_recap_NS$Release_year,
                                                           ncount_recap_NS$recaptureArea),
                                                           '\\.', ' / '), ordered=TRUE)
ncount_recap_NS <- transform(ncount_recap_NS,
                             recaptureArea = factor(recaptureArea, levels = c("W","E","NS")),
                             Tag_area = factor(Tag_area, levels = c(">=54°N","<54°N")))

(p4 <- ggplot(ncount_recap_NS, aes(x =yeararea, fill = Tag_area, y = freq)) +
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label=freq, y = totFreq))+
    #facet_wrap( ~ Release_year, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(face = "bold"),
          #panel.margin = grid::unit(-1.25, "lines"),
          #strip.placement = "outside",
          #strip.background = element_rect(fill = "white", color = "white"),
          axis.text.x = element_text(angle = 90),
          strip.text = element_blank())+
    scale_x_discrete(expand = c(0,0), name = "Release year/area of recapture")+
    scale_y_continuous(expand = c(0,2), "Number of recaptured fish") +
    scale_fill_discrete("Latitude\nRelease"))
(p4 <- ggplot(ncount_recap_NS, aes(x = recaptureArea, fill = Tag_area, y = freq)) +
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label=freq, y = totFreq))+
    facet_wrap( ~ Release_year, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(face = "bold"),
          panel.margin = grid::unit(0, "lines"),
          strip.placement = "outside",
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_blank(),
          axis.text.y = element_blank(),#element_text(angle = 90, hjust = .5),
          strip.background = element_rect(fill = "white", color = "white"))+
    scale_x_discrete(expand = c(.2,.2), name = "Release year")+
    scale_y_continuous(expand = c(0,0), "Number of recaptured fish", breaks = seq(0,200,25)) +
    scale_fill_discrete("Julian\n release"))
# --------------------------

rates_julian <- full_join(ncount_recap_julian %>% rename(freq.recap = freq),
                      ncount_release_julian %>% rename(freq.release = freq) %>%
                        mutate(Release_year = as.integer(as.character(Release_year))),
                      by = c("Release_year", "timing")) %>%
  mutate(rate = freq.recap/freq.release)
rates_julian <- rates_julian %>% group_by(Release_year, recaptureArea) %>% mutate(srate = rate/sum(rate))
names(ncount_recap_julian)
names(ncount_release_julian)
rates_julian <- transform(rates_julian, recaptureArea = factor(recaptureArea, levels = c("W", "E","NS")),
                          timing = factor(timing, levels = c(">=141","<141")))
(p5 <- ggplot(rates_julian, aes(x = recaptureArea, fill = timing, y = srate)) +
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label=freq, y = totFreq))+
    facet_wrap( ~ Release_year, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(face = "bold"),
          panel.margin = grid::unit(0, "lines"),
          strip.placement = "outside",
          #axis.text.x = element_text(size = 12),
          #axis.text.y = element_text(size = 11),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          legend.position = "none",
          axis.text.y = element_text(angle = 90, hjust = .6, size = 11),
          strip.background = element_rect(fill = "white", color = "white"),
          axis.text.x = element_text(angle = 0, size = 12))+
    geom_hline(yintercept = 0.5, lty = 2, lwd = .7)+
    scale_x_discrete(expand = c(.2,.2), name = "Release year",
    )+
    scale_y_continuous(expand = c(0,0.0002), "Recapture rate (% of total rate)",
                       label = scales::percent) +
    scale_fill_discrete("Tagging area"))
#-----
rates_NS <- full_join(ncount_recap_NS %>% rename(freq.recap = freq),
                      ncount_release_NS %>% rename(freq.release = freq) %>%
                        mutate(Release_year = as.integer(as.character(Release_year)),
                               totFreq = NULL),
                      by = c("Release_year", "Tag_area")) %>%
  mutate(rate = freq.recap/freq.release)
rates_NS <- rates_NS %>% group_by(Release_year, recaptureArea) %>% mutate(srate = rate/sum(rate))

names(ncount_recap_NS)
names(ncount_release_NS)
rates_NS <- transform(rates_NS, recaptureArea = factor(recaptureArea, levels = c("W","E","NS")),
                      Tag_area = factor(Tag_area, levels = c(">=54°N","<54°N")))
(p6 <- ggplot(rates_NS, aes(x = recaptureArea, fill = Tag_area, y = srate)) +
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label=freq, y = totFreq))+
    facet_wrap( ~ Release_year, strip.position = "bottom", scales = "free_x", nrow = 1) +
    scale_x_discrete(expand = c(.2,.2), name = "Release year",
    )+
    scale_y_continuous(expand = c(0,0.0002), "Relative recapture rate (%)",
                       label = scales::percent) +
    scale_fill_discrete("Julian\n release")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          #legend.title = element_text(face = "bold"),
          panel.margin = grid::unit(0, "lines"),
          strip.placement = "outside",
          legend.position = "none",
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.y = element_text(angle = 90, hjust = .5),
          strip.background = element_rect(fill = "white", color = "white"),
          axis.text.x = element_text(angle = 0))+
geom_hline(yintercept = 0.5, lty = 2, lwd = .7)
)
# -----------------
# Combination !
# ---------------

library(ggpubr)
g4 <- ggarrange(p1,p2,p3,p4, p5,p6, ncol = 2, nrow = 3)
ggsave(g4, filename = "plots/Arils_supplementary_multipanel_figure_relative_2.tiff",
       width=35, height=30, units="cm", device = "tiff", dpi = "retina")
