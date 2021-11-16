#########################################################################
###
### Content:
### This is the main file to perform some analyses of the mackerel tagging data
###
### Author:
### Kotaro Ono
### Contributor:
###
### Version: 1
###
################################################

#### Set the library path
.libPaths(c("C:/Program Files/R/R-3.6.2/library", "C:/Program Files/R/R-4.0.5/library"))

#### Installing required packages
#devtools::install_github("IMRpelagic/taggart", dependencies = FALSE)
# devtools::install_github("hafro/geo")
# source("R/download_data_functions.R")
source("R/Functions.R")
library(taggart)
library(adehabitatHR)
library(tidyverse)
library(ggmap)
library(raster)
library(mgcv)
library(DHARMa)
library(bbmle) # package to perform model selection using AIC, AICc, etc.
library(glmmTMB)
library(gridExtra)
library(scatterpie)
library(ggplot2)
library(ggnewscale)
library(ggthemes)
library(TMB)
library(TMBhelper)
library(ggpubr) # very nice option for arrange and labeling multipanel plot for manuscript
library(sf)

Norway <- st_read("D:/Dropbox/IMR_projects/Shapefiles/ne_10m_land.shp")
Norway <- st_crop(Norway, c(xmin = -35, ymin = 51, xmax = 21, ymax = 73))


#### Downloading data and checking it
#	tg_catches()         %>% glimpse()
#	tg_catches_bio()     %>% glimpse()
#	tg_expeditions()     %>% glimpse()
#	tg_expeditions_bio() %>% glimpse()
# species <- "mackerel"
# Link_catch_sample <- jsonlite::fromJSON(paste0("http://smartfishsvc.hi.no/api/data/BioSamplesCatches/", species[1]))

Catch_data <- tg_catches()
Catch_bio <- tg_catches_bio()
Mark_recap_data <- tg_expeditions()
Mark_recap_bio <- tg_expeditions_bio()

Data_mackerel <- Mark_recap_data %>% subset(., species == "mackerel")

Catch_data_test <- Catch_data %>% mutate(date = as.Date(ProcessingDate))
Catch_data_test$Year <- as.numeric(substr(Catch_data_test$date, start=1, stop=4))
Catch_data_test <- Catch_data_test %>% filter(Year>=2014)
Catch_data_test$julian_catch <-  as.numeric(julian(as.POSIXct(Catch_data_test$date), as.POSIXct(paste0(2014, "-01-01"), tz = "GMT")))
Catch_data_test$julian_catch_std <-  Catch_data_test$julian_catch %% 365

world <- st_read("shapefile/ne_10m_land.shp")
Atlantic <- st_crop(world, c(xmin = -65, ymin = 46, xmax = 32, ymax = 80))

sum(is.na(Catch_data$cLon))
sum(is.na(Catch_data$cLat))
#map_area <- get_map(location= c(left = -35, bottom = 51, right = 21, top = 72))
#ggmap(map_area) + geom_point(data=Data_mackerel, aes(x=Longitude, y=Latitude))

Catch_data[which(Catch_data$pkid == "2a57fb1b-300e-48ee-b27d-03e4fe1c99fd"),]
Catch_bio[which(Catch_bio$id == "ee6dd054-9274-4ee3-b564-00f233cc941e"),]

#### Define the data frame: creating and merging variables
### Some base variables
Data_mackerel_all <- Data_mackerel %>% left_join(Catch_data %>% dplyr::select(-c('ICES_Rectangle','species','FactoryID','BioSampleID','BioSample')), by=c("CatchID" = "CatchID"))
Data_mackerel_all$dist <- pointDistance(matrix(cbind(Data_mackerel_all$Longitude, Data_mackerel_all$Latitude),ncol=2), matrix(cbind(Data_mackerel_all$cLon, Data_mackerel_all$cLat), ncol=2), lonlat=TRUE)
Data_mackerel_all$duration <- Data_mackerel_all$RecaptureDate - Data_mackerel_all$ReleaseDate
Data_mackerel_all$Nb_years <- as.numeric(lubridate::year(as.Date(Data_mackerel_all$RecaptureDate))) - (lubridate::year(as.Date(Data_mackerel_all$ReleaseDate))) + 1
Data_mackerel_all$Tagger <- as.factor(Data_mackerel_all$Tagger )
Data_mackerel_all$duration <- as.numeric(Data_mackerel_all$duration )
Data_mackerel_all$duration_std <- Data_mackerel_all$duration/Data_mackerel_all$Nb_years
Data_mackerel_all$duration_std1 <- Data_mackerel_all$duration%%365
Data_mackerel_all$Length <- as.numeric(Data_mackerel_all$Length )
Data_mackerel_all$dist <- as.numeric(Data_mackerel_all$dist )
Data_mackerel_all$Release_year <- as.factor(lubridate::year(as.Date(Data_mackerel_all$ReleaseDate)))
Data_mackerel_all$Release_month <- as.factor(lubridate::month(as.Date(Data_mackerel_all$ReleaseDate)))
Data_mackerel_all$Release_month <- factor(Data_mackerel_all$Release_month, labels=c("05","06","08","09","10"))
Data_mackerel_all$pos <- numFactor(Data_mackerel_all$Longitude, Data_mackerel_all$Latitude)
Data_mackerel_all$Catch_year <- as.factor(lubridate::year(as.Date(Data_mackerel_all$RecaptureDate)))
Data_mackerel_all$Catch_month <- as.factor(lubridate::month(as.Date(Data_mackerel_all$RecaptureDate)))
Data_mackerel_all$Release_day <- as.factor(lubridate::day(as.Date(Data_mackerel_all$ReleaseDate)))
Data_mackerel_all$Release_day <- factor(Data_mackerel_all$Release_day, labels=formatC(sort(unique(Data_mackerel_all$Release_day)), width=2, flag=0))
Data_mackerel_all$Release_monthday <- paste(Data_mackerel_all$Release_month, Data_mackerel_all$Release_day, sep="_")
Data_mackerel_all$Release_monthday <- ordered(Data_mackerel_all$Release_monthday)
Data_mackerel_all$ID <- factor(rep(1, nrow(Data_mackerel_all)))
Data_mackerel_all$firstyear <- as.Date(paste0(Data_mackerel_all$Release_year,"-01-01"))
Data_mackerel_all$date <- as.numeric(as.Date(Data_mackerel_all$ReleaseDate) - Data_mackerel_all$firstyear)
Data_mackerel_all$julian_release <-  as.numeric(julian(Data_mackerel_all$ReleaseDate, as.POSIXct(paste0(2014, "-01-01"), tz = "GMT")))
Data_mackerel_all$julian_release_std <-  Data_mackerel_all$julian_release %% 365
Data_mackerel_all$julian_recapture <-  as.numeric(julian(Data_mackerel_all$RecaptureDate, as.POSIXct(paste0(2014, "-06-01"), tz = "GMT")))
Data_mackerel_all$julian_recapture_std <-  Data_mackerel_all$julian_recapture %% 365

## Mvt rate = response variable in the analysis below
Data_mackerel_all$rate_annual <- Data_mackerel_all$dist / Data_mackerel_all$duration_std
Data_mackerel_all$log_rate_annual <- log(Data_mackerel_all$dist/Data_mackerel_all$duration_std)
Data_mackerel_all$rate_annual1 <- Data_mackerel_all$dist / Data_mackerel_all$duration_std1
Data_mackerel_all$log_rate_annual1 <- log(Data_mackerel_all$dist/Data_mackerel_all$duration_std1)
Data_mackerel_all$rate <- Data_mackerel_all$dist / Data_mackerel_all$duration
Data_mackerel_all$log_rate <- log(Data_mackerel_all$dist/Data_mackerel_all$duration)

### Mvt N-S or E-W
Data_mackerel_all$EW_move <- as.numeric(pointDistance(matrix(cbind(Data_mackerel_all$Longitude, Data_mackerel_all$cLat),ncol=2), matrix(cbind(Data_mackerel_all$cLon, Data_mackerel_all$cLat), ncol=2), lonlat=TRUE))
Data_mackerel_all$EW_move <- ifelse(Data_mackerel_all$Longitude < Data_mackerel_all$cLon, Data_mackerel_all$EW_move, -Data_mackerel_all$EW_move)
Data_mackerel_all$NS_move <- as.numeric(pointDistance(matrix(cbind(Data_mackerel_all$Longitude, Data_mackerel_all$Latitude),ncol=2), matrix(cbind(Data_mackerel_all$Longitude, Data_mackerel_all$cLat), ncol=2), lonlat=TRUE))
Data_mackerel_all$NS_move <- ifelse(Data_mackerel_all$Latitude < Data_mackerel_all$cLat, Data_mackerel_all$NS_move, -Data_mackerel_all$NS_move)

### Length bin
Data_mackerel_all$length_bin <- cut(Data_mackerel_all$Length, breaks=c(0, 34, 37, 45))
Data_mackerel_all$length_bin  <- factor(Data_mackerel_all$length_bin, labels=c("(0,34cm]", "(34cm,37cm]", "(37cm,45cm]"))

### Need to make some geographical division of the data
South_Ireland <- unique(Data_mackerel_all$ICES_Rectangle)[unlist(sapply(31:36, function(x) grep(x, unique(Data_mackerel_all$ICES_Rectangle))))]
North_Ireland <- unique(Data_mackerel_all$ICES_Rectangle)[unlist(sapply(37:46, function(x) grep(x, unique(Data_mackerel_all$ICES_Rectangle))))]
Bergen <- unique(Data_mackerel_all$ICES_Rectangle)[grep("F", unique(Data_mackerel_all$ICES_Rectangle))]
Iceland <- unique(Data_mackerel_all$ICES_Rectangle)[grep(58, unique(Data_mackerel_all$ICES_Rectangle))]
Data_mackerel_all <- Data_mackerel_all %>% mutate(Tag_area = ifelse(ICES_Rectangle %in% South_Ireland, "South_Ireland",
                                                                    ifelse(ICES_Rectangle %in% North_Ireland, "North_Ireland",
                                                                           ifelse(ICES_Rectangle %in% Bergen, "Bergen",
                                                                                  ifelse(ICES_Rectangle %in% Iceland, "Iceland", "Weird")))))
Data_mackerel_all <- Data_mackerel_all %>% mutate(Tag_area_large = ifelse(ICES_Rectangle %in% c(South_Ireland,North_Ireland), "Ireland",
                                                                          ifelse(ICES_Rectangle %in% Bergen, "Bergen",
                                                                                 ifelse(ICES_Rectangle %in% Iceland, "Iceland", "Weird"))))

Data_mackerel_all$Tag_area <- as.factor(Data_mackerel_all$Tag_area)
Data_mackerel_all$Tag_area_large <- as.factor(Data_mackerel_all$Tag_area_large)

Data_mackerel_all$Release_timing <- ifelse(Data_mackerel_all$Release_monthday < "05_22", "First_half", "Second_half")
#Data_mackerel_all$Release_timing[which(Data_mackerel_all$Release_monthday >= "05_16" & Data_mackerel_all$Release_monthday <= "05_26")] <- "Mid"

Data_mackerel_all$Category <- paste(Data_mackerel_all$Tag_area, Data_mackerel_all$Release_timing, sep="_")

#### If I want to look at both the mark and recapture
Data_mackerel_final <- Data_mackerel_all %>% subset(., !is.na(dist))

Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year) %>%
  mutate(mean_lon_tag=mean(Longitude),
         mean_lat_tag=mean(Latitude),
         mean_lon_recap=mean(cLon),
         mean_lat_recap=mean(cLat))
Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year, length_bin) %>%
  mutate(mean_lon_recap_size=mean(cLon),
         mean_lat_recap_size=mean(cLat))
Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year, length_bin, Tag_area_large) %>%
  mutate(mean_lon_recap_size_area=mean(cLon),
         mean_lat_recap_size_area=mean(cLat))
Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year, Tag_area_large) %>%
  mutate(mean_lon_tag_area=mean(Longitude),
         mean_lat_tag_area=mean(Latitude))

Data_mackerel_use_Ireland <- subset(Data_mackerel_final, Tag_area %in% c("South_Ireland", "North_Ireland"))

## Using duration (not standardized) but focusing on within year recapture
# Data_mackerel_use_Ireland_select <- subset(Data_mackerel_use_Ireland, duration < 180 & Release_year%in%2014:2020)
# Data_mackerel_use_Ireland_select$Release_year <- as.factor(as.character(Data_mackerel_use_Ireland_select$Release_year))

## Do some explanatory analysis to explore what is causing these patterns of residuals
mean_tag_area <- Data_mackerel_use_Ireland_select %>% group_by(Tag_area) %>% summarize(median = median(log_rate))
mean_diff_tag_area <- as.numeric(abs(mean_tag_area[1,2] - mean_tag_area[2,2]))

Skip_gam_analysis = TRUE
do_plotting = FALSE
do_bayesian = FALSE
model_choice = "lme"   # choice of "lm" or "lme"


if (do_plotting == TRUE){

  #### Some visual inspection of the tagging data
  ### How often do mackerel move south
  ggplot(Data_mackerel_all[which((Data_mackerel_all$Latitude > Data_mackerel_all$cLat) ==TRUE),]) +
    geom_point(aes(x=Longitude, y=Latitude)) + geom_point(aes(x=cLon, y=cLat), col="red")

  ### Something is going on.... let's dig a bit more
  ## Release month and location
  with(Data_mackerel_final, table(Release_month, Tag_area_large))
  with(Data_mackerel_final, table(Release_month, Tag_area_large, Release_year))
  with(Data_mackerel_final, table(Release_month, Tag_area))
  with(Data_mackerel_final, table(Release_month, Tag_area, Release_year))

  ## Catch month and location
  with(Data_mackerel_final, table(Catch_month, ICES_Rectangle))
  with(Data_mackerel_final, table(Catch_month, recatch_ICES_Rectangle, Catch_year))
  with(Data_mackerel_final, table(Catch_month, Tag_area))
  with(Data_mackerel_final, table(Catch_month, Tag_area, Release_year))

  ## Release month and length
  ggplot(Data_mackerel_final, aes(x=Release_month, y=length), alpha=0.5) + geom_boxplot() +
    facet_grid(Tag_area_large~ Release_year) + theme_bw() #+ coord_flip()

  ## the duration
  ggplot(Data_mackerel_final, aes(x=length_bin, y=duration)) + facet_grid(Tag_area_large ~ Release_year) +
    geom_boxplot() + theme_bw()

  ggplot(Data_mackerel_final, aes(x=length_bin, y=duration_std)) + facet_grid(Tag_area_large ~ Release_year) +
    geom_boxplot() + theme_bw()

  ggplot(Data_mackerel_final, aes(x=length_bin, y=duration_std1)) + facet_grid(Tag_area_large ~ Release_year) +
    geom_boxplot() + theme_bw()

  # travel distance
  (ncount_year <- Data_mackerel_final %>% group_by(Release_year,length_bin) %>% count())
  ggplot(Data_mackerel_final, aes(x=length_bin, y=dist/1000)) + facet_grid(Tag_area_large ~ Release_year) +
    geom_boxplot() + theme_bw()

  # travel rate (linear but still)
  ncount_year <- Data_mackerel_final %>% group_by(Release_year) %>% count()
  ggplot(Data_mackerel_final, aes(x=length_bin, y=log_rate_annual1)) + facet_grid(Tag_area_large ~ Release_year) +
    geom_boxplot() + theme_bw()

  # length
  with(Data_mackerel_final, boxplot(length ~ Release_year, ylim=c(0,50), ylab="Mackerel size (cm)"))
  ncount_year <- Data_mackerel_final %>% group_by(Release_year) %>% count()
  text(1:9, rep(48, 9), sapply(1:9, function(x) paste0("n=", ncount_year[x,2])), cex=0.9)


  Data_mackerel_use_Ireland_select
  with(Data_mackerel_use_Ireland_select, boxplot(length ~ Release_year, ylim=c(0,50), ylab="Mackerel size (cm)"))
  ncount_year <- Data_mackerel_use_Ireland_select %>% group_by(Release_year) %>% count()
  text(1:9, rep(48, 9), sapply(1:9, function(x) paste0("n=", ncount_year[x,2])), cex=0.9)

  #### Some plotting of the tagging location and recapture
  #### Maps + centre of gravity
  ### By Year
  gg <- list()
  for (ijk in seq_along(2011:2019))
  {
    gg[[ijk]] <- ggmap(map_area) + geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=Longitude, y=Latitude)) +
      + ggtitle((2011:2019)[ijk]) + theme(plot.title = element_text(hjust = 0.5)) +
      geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=cLon, y=cLat, col=length_bin), alpha=0.5) +
      geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_recap_size, y=mean_lat_recap_size, col=length_bin), shape=7, size=5) +
      geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_tag, y=mean_lat_tag), col="black", shape=7, size=5)

  }

  grid.arrange(gg[[1]],gg[[2]],gg[[3]],
               gg[[4]],gg[[5]],gg[[6]],
               gg[[7]],gg[[8]],gg[[9]],
               nrow=3, ncol=3)

  ### by year, region and size
  ## all years
  ggg <- list()
  for (ijk in seq_along(2011:2019))
  {
    ggg[[ijk]] <- ggmap(map_area) + geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=Longitude, y=Latitude)) +
      facet_grid(. ~ Tag_area_large) + ggtitle((2011:2019)[ijk]) + theme(plot.title = element_text(hjust = 0.5)) +
      geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=cLon, y=cLat, col=length_bin), alpha=0.5) +
      geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_recap_size_area, y=mean_lat_recap_size_area, col=length_bin), shape=7, size=5) +
      geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_tag_area, y=mean_lat_tag_area), col="black", shape=7, size=5)
  }
  grid.arrange(ggg[[1]],ggg[[2]],ggg[[3]],
               nrow=3, ncol=1)
  grid.arrange(ggg[[4]],ggg[[5]],ggg[[6]],
               nrow=3, ncol=1)
  grid.arrange(ggg[[7]],ggg[[8]],ggg[[9]],
               nrow=3, ncol=1)

  ## only 1 year subset
  gggg <- list()
  for (ijk in seq_along(2011:2019))
  {
    gggg[[ijk]] <- ggmap(map_area) + geom_point(data=subset(Data_mackerel_final_sameyear, Release_year==(2011:2019)[ijk]), aes(x=Longitude, y=Latitude)) +
      facet_grid(. ~ Tag_area_large) + ggtitle((2011:2019)[ijk]) + theme(plot.title = element_text(hjust = 0.5)) +
      geom_point(data=subset(Data_mackerel_final_sameyear, Release_year==(2011:2019)[ijk]), aes(x=cLon, y=cLat, col=length_bin), alpha=0.5) +
      geom_point(data=subset(Data_mackerel_final_sameyear, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_recap_size_area, y=mean_lat_recap_size_area, col=length_bin), shape=7, size=5) +
      geom_point(data=subset(Data_mackerel_final_sameyear, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_tag_area, y=mean_lat_tag_area), col="black", shape=7, size=5)
  }
  grid.arrange(gggg[[1]],gggg[[2]],gggg[[3]],
               nrow=3, ncol=1)
  grid.arrange(gggg[[4]],gggg[[5]],gggg[[6]],
               nrow=3, ncol=1)
  grid.arrange(gggg[[7]],gggg[[8]],gggg[[9]],
               nrow=3, ncol=1)



}

########## Preparing the final data from Ireland for the analysis

Data_mackerel_use_Ireland1 <- Data_mackerel_use_Ireland
Data_mackerel_use_Ireland1$Catch_month <- as.numeric(as.character(Data_mackerel_use_Ireland1$Catch_month))
Data_mackerel_use_Ireland1$Catch_year <- as.numeric(as.character(Data_mackerel_use_Ireland1$Catch_year))
Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(1,2)),'Catch_month'] <- Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(1,2)),'Catch_month'] + 12
Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(13,14)),'Catch_year'] <- Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(13,14)),'Catch_year'] - 1
Data_mackerel_use_Ireland_select_origin <- Data_mackerel_use_Ireland1 %>%
  filter(Catch_month %in% c(7:14), as.numeric(as.character(Catch_year)) == as.numeric(as.character(Release_year)), Release_year%in%2014:2020)
Data_mackerel_use_Ireland_select_origin_year1 <- Data_mackerel_use_Ireland1 %>%
  filter(Catch_month %in% c(7:14), as.numeric(as.character(Catch_year)) == as.numeric(as.character(Release_year))+1, Release_year%in%2014:2020)
Data_mackerel_use_Ireland_select_origin_year2 <- Data_mackerel_use_Ireland1 %>%
  filter(Catch_month %in% c(7:14),  as.numeric(as.character(Catch_year)) == as.numeric(as.character(Release_year))+2, Release_year%in%2014:2020)
Data_mackerel_use_Ireland_select_origin$Catch_year <- as.factor(Data_mackerel_use_Ireland_select_origin$Catch_year)
Data_mackerel_use_Ireland_select_origin_year1$Catch_year <- as.factor(Data_mackerel_use_Ireland_select_origin_year1$Catch_year)
Data_mackerel_use_Ireland_select_origin_year2$Catch_year <- as.factor(Data_mackerel_use_Ireland_select_origin_year2$Catch_year)


ggplot(Norway) + geom_sf() + geom_jitter(data=Data_mackerel_use_Ireland_select_origin, aes(x=cLon, y=cLat, col=factor(Catch_month))) +
  theme_bw() + facet_wrap(~Catch_year) + geom_hline(yintercept=62) + geom_vline(xintercept=-10)+ geom_vline(xintercept=-4, col="red")


month = 12
best="overall"
data_origin="0year"

run_directionality <- function(month, best="overall", data_origin="0year", parallel=FALSE){
  ## Doing the analysis of P(moving to Iceland) or P(moving to Norway)
  Limit_month <- month ## c(9,10,11)
  if (month == 9)   label <- "Sep30th"
  if (month == 10)  label <- "Oct31st"
  if (month == 11)  label <- "Nov30st"
  if (month == 12)  label <- "Dec31st"
  if (month == 13)  label <- "Jan31st"
  if (month == 14)  label <- "Feb28th"

  # Full data
  if (data_origin=="0year") data = Data_mackerel_use_Ireland_select_origin
  if (data_origin=="1year") data = Data_mackerel_use_Ireland_select_origin_year1
  if (data_origin=="2year") data = Data_mackerel_use_Ireland_select_origin_year2
  data$month <- as.numeric(data$Catch_month)
  Data_mackerel_use_Ireland_select <- subset(data, Release_year %in% 2014:2020) %>% filter(month <= Limit_month)
  # Selected data
  # Data_mackerel_use_Ireland_select_filter <- subset(Data_mackerel_use_Ireland_select, cLon>-10)
  #Data_mackerel_use_Ireland_select_filter <- subset(Data_mackerel_use_Ireland_select, cLon< -10)

  Data_mackerel_use_Ireland_select$ID <- 1:nrow(Data_mackerel_use_Ireland_select)

  # Rescaling parameters to ease interpretation
  Data_mackerel_use_Ireland_select$julian_recapture_scaled <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std)
  Data_mackerel_use_Ireland_select$julian_recapture_standardized <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std, center=FALSE)
  Data_mackerel_use_Ireland_select$length_scaled <- scale(Data_mackerel_use_Ireland_select$Length)
  Data_mackerel_use_Ireland_select$Latitude_scaled <- scale(Data_mackerel_use_Ireland_select$Latitude)
  Data_mackerel_use_Ireland_select$Latitude2_scaled <- scale((Data_mackerel_use_Ireland_select$Latitude)^2)

  Data_mackerel_use_Ireland_select <- Data_mackerel_use_Ireland_select %>%
    mutate(toiceland= ifelse(cLon < -10 & cLat>=62, 1, 0),
           dist_toiceland = abs(as.numeric(pointDistance(matrix(cbind(max(cLon, na.rm=T)+0.01, cLat),ncol=2), matrix(cbind(cLon, cLat), ncol=2), lonlat=TRUE))),
           to_norway = ifelse((cLon >= -10 & cLat>=62), 1, 0),
           to_northsea = ifelse((cLon >= -10 & cLat<62), 1, 0),
           direction = ifelse(cLon < -10 & cLat>=62, 1L, ifelse((cLon >= -10 & cLat>=62), 2L, 0L)))

  Data_mackerel_use_Ireland_select %>% group_by(toiceland, Catch_year) %>% summarize(n=n())

  Data_mackerel_use_Ireland_select %>% filter(toiceland == 1) %>% ggplot(aes(x=Latitude)) + geom_histogram()


  hist(Data_mackerel_use_Ireland_select$dist_toiceland, breaks=50, xlab="Distance from Norway coast (m)", main="")

  ## Analysis of the directional movement (multinomial)
  www <- Data_mackerel_use_Ireland_select
  www <- www %>% mutate(group = paste0(Catch_year, "_", Release_timing),
                        group = as.factor(group),
                        Release_timing_fact = as.factor(Release_timing),
                        group = droplevels(group),
                        Catch_year = droplevels(Catch_year),
                        julian_recapture_std_scaled = scale(julian_recapture_std),
                        julian_recapture_scaled = scale(julian_recapture),
                        Latitude_scaled = scale(Latitude),
                        Length_scaled = scale(Length),
                        toiceland_bin = as.factor(toiceland),
                        tonorway_bin = as.factor(to_norway),
                        tonorthsea_bin = as.factor(to_northsea)
  )
  m0 <- gam(list(direction ~  1, ~  1),family=multinom(K=2), data = www)
  # m00 <- gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) ,
  #                 ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled)) ,family=multinom(K=2), data = www)
  # m01 <- gam(list(direction ~  Catch_year+s(julian_recapture_std_scaled), ~  Catch_year+s(julian_recapture_std_scaled)),family=multinom(K=2), data = www)
  # m02 <- gam(list(direction ~  Catch_year+ s(Length_scaled)+s(julian_recapture_std_scaled), ~  Catch_year+ s(Length_scaled)+s(julian_recapture_std_scaled)),family=multinom(K=2), data = www)
  m1 <- gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                 ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled)),
            family=multinom(K=2), data = www)
  m2 <- gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                 ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled, by=Catch_year) + s(julian_recapture_std_scaled)),
            family=multinom(K=2), data = www)
  m3 <- gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                 ~  -1 + s(Latitude_scaled, by=Catch_year) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled)),
            family=multinom(K=2), data = www)
  m4 <- gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                 ~  -1 + s(Latitude_scaled, by=Catch_year) + Catch_year + s(Length_scaled, by=Catch_year) + s(julian_recapture_std_scaled)),
            family=multinom(K=2), data = www)
  # m4 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled, by=Catch_year)),
  #           family=multinom(K=2), data = www)
  # m5 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled, by=Catch_year)),
  #           family=multinom(K=2), data = www)
  # m6 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled, by=Catch_year)),
  #           family=multinom(K=2), data = www)
  # m8 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled, by=Catch_year)),
  #           family=multinom(K=2), data = www)
  # m9 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled, by=Catch_year)),
  #           family=multinom(K=2), data = www)
  # m3 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude) + Catch_year +  s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www) # failed
  # m4 <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # m5 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # m6 <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # m7 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # m8 <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # mall <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # mall1 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  # mall2 <- gam(list(direction ~  -1 + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  #                ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  #           family=multinom(K=2), data = www)
  (AICs <- AIC(m0,m1,m2,m3,m4))
  (BICs <- BIC(m0,m1,m2,m3,m4))

  # plot.gam(m1, pch=19, cex=0.5, scale=0,pages = 1)
  # library(generalhoslem)
  # pred <- predict(m2, type="response")
  # logitgof(as.factor(www$direction), pred)
  # pred <- predict(m6, type="response")
  # logitgof(as.factor(www$direction), pred)
  # pred <- predict(m1, type="response")
  # logitgof(as.factor(www$direction), pred)
  #


  ### 10-fold cross-validation

  library(groupdata2)
  library(furrr)
  library(progressr)

  set.seed(777)
  kfolding <- function(it){
     p()
     ERRORS <- c()
     www1 <- www %>% ungroup() %>% fold(k = 10)
      for (k in 1:10){
        testing <- www1[www1$.folds == k,]
        training <- www1[-testing$ID,]

        m0=m1=m2=m3=m4=NULL
        m0 <- gam(list(direction ~  1, ~  1),family=multinom(K=2), data = training)
        m1 <- try(gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                           ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled)),
                      family=multinom(K=2), data = training))
        m2 <- try(gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                           ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled, by=Catch_year) + s(julian_recapture_std_scaled)),
                      family=multinom(K=2), data = training))
        m3 <- try(gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                           ~  -1 + s(Latitude_scaled, by=Catch_year) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled)),
                      family=multinom(K=2), data = training))
        m4 <- try(gam(list(direction ~  -1 + s(Latitude_scaled) + Catch_year + s(Length_scaled) + s(julian_recapture_std_scaled),
                           ~  -1 + s(Latitude_scaled, by=Catch_year) + Catch_year + s(Length_scaled, by=Catch_year) + s(julian_recapture_std_scaled)),
                      family=multinom(K=2), data = training))

        mse0 = mse1 = mse2 = mse3 = mse4 = NA

        if (is(m0) != "try-error")   pred0 <- predict(m0, type="response", newdata = testing)
        if (is(m1) != "try-error")   pred1 <- predict(m1, type="response", newdata = testing)
        if (is(m2) != "try-error")   pred2 <- predict(m2, type="response", newdata = testing)
        if (is(m3) != "try-error")   pred3 <- predict(m3, type="response", newdata = testing)
        if (is(m4) != "try-error")   pred4 <- predict(m4, type="response", newdata = testing)

        tonum <- function(x) as.numeric(as.character(x))
        obs <- data.frame(other=tonum(testing$tonorthsea_bin), Ice=tonum(testing$toiceland_bin), Nor=tonum(testing$tonorway_bin))
        if (is(m0) != "try-error")   mse0 <- mean(apply(pred0-obs, 1, function(x) sum(x^2)))
        if (is(m1) != "try-error")   mse1 <- mean(apply(pred1-obs, 1, function(x) sum(x^2)))
        if (is(m2) != "try-error")   mse2 <- mean(apply(pred2-obs, 1, function(x) sum(x^2)))
        if (is(m3) != "try-error")   mse3 <- mean(apply(pred3-obs, 1, function(x) sum(x^2)))
        if (is(m4) != "try-error")   mse4 <- mean(apply(pred4-obs, 1, function(x) sum(x^2)))

        error <- c(mse0, mse1, mse2, mse3, mse4)
        ERRORS <- rbind(ERRORS, error)
      }
     return(ERRORS)
    }

  if (parallel == TRUE){
    future::plan("multicore")
    system.time(
      with_progress({
      p <- progressor(steps = 30)
      ERRORS <- furrr::future_pmap_dfr(data.frame(it=1:30), kfolding,
                                    .options = furrr::furrr_options(seed = TRUE))
    })
    )
  }

  if (parallel == FALSE){
    ERRORS <- c()
    for (it in 1:30) {
      err <- kfolding(it)
      ERRORS <- rbind(ERRORS, err)
    }
    print(it)
  }

  Res_cross_validation <- apply(ERRORS,2, mean, na.rm=T)

  Ranking <- data.frame(AIC=rank(AICs$AIC), BIC=rank(BICs$BIC), CV=rank(Res_cross_validation))

  if (best == "overall") {
    which.best = which.min(apply(Ranking, 1, sum))
    best_model <- list(m0,m1,m2,m3,m4)[which.best]
  }


  out <- list()
  out$scenario <- paste0("month=", month, "_lag=", data_origin)
  out$best_model_id <- which.best
  out$best_model <- best_model

  saveRDS(out, file=paste0("month=", month, "_lag=", data_origin, ".rds"))


  return(out)

}


Dec_0lag <- run_directionality(month=12, data_origin="0year")
Dec_1lag <- run_directionality(month=12, data_origin="1year")
Dec_2lag <- run_directionality(month=12, data_origin="2year")
Nov_0lag <- run_directionality(month=11, data_origin="0year")
Nov_1lag <- run_directionality(month=11, data_origin="1year")
Nov_2lag <- run_directionality(month=11, data_origin="2year")
Jan_0lag <- run_directionality(month=13, data_origin="0year")
Jan_1lag <- run_directionality(month=13, data_origin="1year")
Jan_2lag <- run_directionality(month=13, data_origin="2year")


  best_model <- m1
  pp <- plot.gam(best_model, pch=19, cex=0.5, scale=0,pages = 1)
  pred <- predict(best_model, type="response")

  # #### marginal effect in real scale
  # datpred <- expand.grid(Latitude_scaled = seq(-2,3,by=0.01), Catch_year=levels(www$Catch_year)[1], Length_scaled=0, julian_recapture_std_scaled=0)
  # pred_marginal <- predict(best_model, newdata=datpred, type="response", se=TRUE)
  # datpred <- data.frame(lat = datpred$Latitude_scaled, y = pred_marginal$fit[,2], se = pred_marginal$se.fit[,2])
  # datpred <- datpred %>% mutate(low = y - 1.96*se, upr=y+1.96*se)
  # ggplot(datpred, aes(x=lat, y=y)) + geom_ribbon(aes(ymin=low, ymax=upr)) + geom_line()
  #
  #
  #
  # datpred$pred1 <- pred_marginal$fit[,1]
  # datpred$pred2 <- pred_marginal$fit[,2]
  # datpred$pred3 <- pred_marginal$fit[,3]
  # datpred$sered1 <- pred_marginal$se.fit[,1]
  # datpred$sepred2 <- pred_marginal$se.fit[,2]
  # datpred$sepred3 <- pred_marginal$se.fit[,3]
  # datpred1 <- reshape2::melt(datpred, measure.vars=5:7, variable.name="direction", id.vars=1:4)
  # datpred2 <- reshape2::melt(datpred, measure.vars=8:10, variable.name="direction", id.vars=1:4, value.name="se")
  # datpred <- datpred1 %>% left_join(datpred2)
  # ggplot(datpred %>% filter(direction != "pred1"), aes(x=Latitude_scaled, y = value)) +
  #   geom_line() + gg_control + facet_grid(direction~.)


  #### Estimated coefficients
  Variable <- c("Iceland: latitude effect", "Iceland: fish body size effect", "Iceland: recapture date effect",
                "Norway: latitude effect", "Norway: fish body size effect", "Norway: recapture date effect")
  Region <- c("Go to Iceland","Go to Iceland","Go to Iceland","Go to Norway","Go to Norway","Go to Norway")
  Variable <- c("Latitude effect","Fish body size effect","Recapture date effect",
                "Latitude effect","Fish body size effect","Recapture date effect")
  Dat_all <- c()
  for (i in 1:6) {
    dat <- data.frame(x=pp[[i]]$x,y=pp[[i]]$fit,se=pp[[i]]$se,ymin=pp[[i]]$fit-1.96*pp[[i]]$se,ymax=pp[[i]]$fit+1.96*pp[[i]]$se)
    dat$Variable <- Variable[i]
    dat$Region <- Region[i]
    Dat_all <- rbind(Dat_all, dat)
  }
  Dat_raw_all <- c()
  for (i in 1:6) {
    dat <- data.frame(x=pp[[i]]$raw, direction=www$direction, y=min(pp[[i]]$fit-1.96*pp[[i]]$se))
    dat$Variable <- Variable[i]
    dat$Region <- Region[i]
    Dat_raw_all <- rbind(Dat_raw_all, dat)
  }


  #### Some plotting of results
  gg_control <- theme(axis.title = element_text(size=15),
                      axis.text = element_text(size=14),
                      plot.title = element_text(hjust=0.5, size=18),
                      strip.text = element_text(size=14))

  p1 <- ggplot(Dat_all, aes(x=x,y=y)) + geom_line(size=2) +
    geom_line(data=Dat_all, aes(x=x, y=ymin), linetype=2) +
    geom_line(data=Dat_all, aes(x=x, y=ymax), linetype=2) +
    facet_wrap(Region~Variable, scales="free") +
    theme_bw() + gg_control + labs(x="Standardized value", y="Estimated effects") +
    geom_rug(data=Dat_raw_all, aes(x=x,y=y))



  library(jtools)
  Catch_year <- data.frame(fit=summary(m1)$p.coeff, se=summary(m1)$se[grep("Catch_year", names(summary(m1)$se))])
  Catch_year$Year <- str_sub(rownames(Catch_year), 11,14)
  Catch_year$Direction <- c(rep("Iceland", nrow(Catch_year)/2), rep("Norway", nrow(Catch_year)/2))
  p2 <- ggplot(Catch_year, aes(x=Year,  y=fit)) + geom_errorbar(aes(ymin=fit-1.96*se, ymax=fit+1.96*se)) + facet_wrap(~Direction, nrow=2, scales="free") +
    theme_bw() + gg_control + labs(y="") + theme(axis.text.x = element_text(size=10))

  ppp <- ggarrange(p1,p2,ncol=2, widths=c(0.8,0.4))
  ggsave(ppp, file=paste0(getwd(), "/MS/figs/Marginal_effects_", month, "_", data_origin, ".pdf"),width=30, height=18, units="cm", dpi = 450)


