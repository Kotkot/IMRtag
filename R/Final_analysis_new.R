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
  tonum <- function(x) as.numeric(as.character(x))

## loadig some shapefiles
  new_proj <- 3035
  Norway <- st_read("D:/Dropbox/IMR_projects/Shapefiles/ne_10m_land.shp")
  Norway <- st_crop(Norway, c(xmin = -35, ymin = 51, xmax = 21, ymax = 73))
  Norway_proj <- st_transform(Norway, crs=new_proj)
  ICESecoregion <- st_read("D:/Dropbox/IMR_projects/Shapefiles/ICES_ecoregions_20171207_erase_ESRI.shp")
  ICESecoregion_cut <- st_crop(ICESecoregion, c(xmin = -35, ymin = 51, xmax = 21, ymax = 73))
  ICESecoregion_proj <- st_transform(ICESecoregion, crs=new_proj)
  ICESecoregion_cut_proj <- st_transform(ICESecoregion_cut, crs=new_proj)

  ggplot(ICESecoregion_cut_proj) + geom_sf(aes(geometry = geometry, fill=Ecoregion))

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

  if (Skip_gam_analysis == FALSE){

      #### Calculate linear distance between mark and recapture location
      #### Then standardize by number of days in between
      #### Model this with respect to the different covariates

      ### Analysis over the whole area: but area South&north of Ireland is the only place with consistent release over the years

      ## Using duration (not standardized)
      # m01 <- gam(log_rate ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_final)
      # m02 <- gam(log_rate ~ s(length) + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_final)
      # m03 <- gam(log_rate ~ s(length) + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_final)
      # m04 <- gam(rate ~ s(length) + s(Longitude, Latitude) + Release_year, family=gaussian(link="log"), data=Data_mackerel_final)
      # m05 <- gam(rate ~ s(length) + s(Longitude, Latitude) + Release_year, family=Gamma(link = "log"), data=Data_mackerel_final)
      # AICtab(m01,m02,m03,m04,m05)
      # plot(m03, all.terms=TRUE, residuals=TRUE, pages=1)
      # summary(m03)

      ## Using duration (not standardized) but focusing only on within year recapture
      m0 <- gam(log_rate ~ s(length) + s(Longitude,Latitude), family=gaussian, data=Data_mackerel_use_select)
      m1 <- gam(log_rate ~ s(length) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_select)
      m2 <- gam(log_rate ~ s(length) + Release_month + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_select)
      m3 <- gam(log_rate ~ s(length) + Release_month + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_select)
      m4 <- gam(log_rate ~ s(length, by= Release_year) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_select)
      m5 <- gam(log_rate ~ s(Longitude, Latitude, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_select)
      m6 <- gam(log_rate ~ s(Longitude, Latitude, by= Tag_area) + Release_year + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_select)
      AICtab(m1,m2,m3,m4,m5,m6)
      plot(m5, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_select <- simulateResiduals(fittedModel = m5, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_select, rank=T, quantreg = TRUE)
      testResiduals(sim_select)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_select$Length, sim_select$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_select$Latitude, sim_select$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_select$Longitude, sim_select$scaledResiduals, main = "Longitude")

      ## Using duration_std
      # mm1 <- gam(log_rate_annual ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_final)
      # mm2 <- gam(log_rate_annual ~ s(length) + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_final)
      # mm3 <- gam(log_rate_annual ~ s(length) + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_final)
      # mm4 <- gam(rate_annual ~ s(length) + s(Longitude, Latitude) + Release_year, family=gaussian(link="log"), data=Data_mackerel_final)
      # mm5 <- gam(rate_annual ~ s(length) + s(Longitude, Latitude) + Release_year, family=Gamma(link = "log"), data=Data_mackerel_final)
      # mm6 <- gam(rate_annual ~ s(length) + s(Longitude, Latitude, by=Tag_area) + Release_year, family=Gamma(link = "log"), data=Data_mackerel_final)
      # mm0 <- gam(log_rate_annual ~ s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_final)
      #
      # AICtab(mm0, mm1,mm2,mm3,mm4,mm5)
      # plot(mm3, all.terms=TRUE, residuals=TRUE, pages=1)

      ## Using duration std1
      # mmm0 <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_final)
      # mmm1 <- gam(log_rate_annual1 ~ s(length) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_final)
      # mmm2 <- gam(log_rate_annual1 ~ s(length) + Release_month + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_final)
      # mmm3 <- gam(log_rate_annual1 ~ s(length) + Release_month + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_final)
      # mmm4 <- gam(log_rate_annual1 ~ s(length, by= Release_year) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_final)
      # mmm5 <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_final)
      # mmm6 <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Tag_area) + s(length) + Release_year + Release_month , family=gaussian, data=Data_mackerel_final)
      # mmm7 <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Tag_area)  + s(length, by= Release_year) + Release_year + Release_month , family=gaussian, data=Data_mackerel_final)
      # mmm7a <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Tag_area)  + s(length, by= Release_year) + Release_month , family=gaussian, data=Data_mackerel_final)
      # mmm7b <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Tag_area)  + s(length, by= Release_year) + Release_year , family=gaussian, data=Data_mackerel_final)
      # mmm7c <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Tag_area)  + s(length, by= Release_year) , family=gaussian, data=Data_mackerel_final)

      # AICtab(mmm0, mmm1,mmm2,mmm3,mmm4,mmm5,mmm6,mmm7,mmm7a,mmm7b,mmm7c)
      # plot(mmm7, all.terms=TRUE, residuals=TRUE)
      # sim_std1 <- simulateResiduals(fittedModel = mmm4, n = 1000, integerResponse = FALSE, plot=FALSE)
      # plot(sim_std1, rank=T, quantreg = TRUE)
      # testResiduals(sim_std1)
      # par(mfrow=c(3,2))
      # plotResiduals(Data_mackerel_final$Length, sim_std1$scaledResiduals, main = "length")
      # plotResiduals(Data_mackerel_final$Release_month, sim_std1$scaledResiduals, main = "Release_month")
      # plotResiduals(Data_mackerel_final$Release_year, sim_std1$scaledResiduals, main = "Release_year")
      # plotResiduals(Data_mackerel_final$Latitude, sim_std1$scaledResiduals, main = "Latitude")
      # plotResiduals(Data_mackerel_final$Longitude, sim_std1$scaledResiduals, main = "Longitude")

      #m1a <- gamm(rate ~ s(length) + s(Longitude, Latitude), random=list(tagger= ~1), family=Gamma, data=Data_mackerel_final)
      #m2a <- gamm(rate ~ s(length) + s(Longitude, Latitude) + s(date), random=list(tagger= ~1), family=Gamma, data=Data_mackerel_final)
      #mm1 <- glmmTMB(log_rate ~ length + mat(pos+0|ID) + Release_year, family=gaussian, data=Data_mackerel_final)


      ### Limiting the analysis to releases around Ireland:


      with(Data_mackerel_use_Ireland_select, plot(length, rate))
      with(Data_mackerel_use_Ireland_select, plot(Tag_area, rate))
      with(Data_mackerel_use_Ireland_select, plot(Release_timing, rate))
      with(Data_mackerel_use_Ireland_select, plot(julian_recapture_std, rate))


      m00 <- gam(log_rate ~ s(length_scaled) + s(Latitude) + factor(Release_year) + julian_recapture_scaled, data=Data_mackerel_use_Ireland_select)
      m01 <- gam(log_rate ~ s(length) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
      m02 <- gam(log_rate ~ s(length) + Category, family=gaussian, data=Data_mackerel_use_Ireland_select)
      m03 <- gam(log_rate ~ s(length) + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m04 <- gam(log_rate ~ s(length) + s(Longitude, Latitude) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m05 <- gam(log_rate ~ length + s(Longitude, Latitude, by= Release_year) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m06 <- gam(log_rate ~ length + s(Longitude, Latitude) + Release_year + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m07 <- gam(log_rate ~ length + s(Longitude, Latitude) + Release_year + s(julian_recapture_std, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)

      m0 <- gam(log_rate ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m1 <- gam(log_rate ~ s(length) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m2 <- gam(log_rate ~ s(length) + Release_month + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m3 <- gam(log_rate ~ s(length) + Release_month + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
      m4 <- gam(log_rate ~ s(length, by= Release_year) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m5 <- gam(log_rate ~ s(Longitude, Latitude, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland_select)
      m6 <- gam(log_rate ~ s(Longitude, Latitude, by= Release_year) + s(length) , family=gaussian, data=Data_mackerel_use_Ireland_select)
      m7 <- gam(log_rate ~ s(Longitude, Latitude, by= Release_year) + length , family=gaussian, data=Data_mackerel_use_Ireland_select)
      m8 <- gam(log_rate ~ s(Longitude, Latitude, by= Release_year) + s(length, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)

      m10 <- gam(log_rate ~ length + s(Longitude, Latitude) + s(date, by=Release_year) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
      m11 <- gam(log_rate ~ length + s(Longitude, Latitude) + s(date, by=Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m12 <- gam(log_rate ~ length + s(Longitude, Latitude, by=Release_year) + s(date, by=Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m13 <- gam(log_rate ~ length + s(Longitude, Latitude, by=Release_year) + s(julian_recapture_std, by=Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)


      library(glmmTMB)
      library(DHARMa)
      m14 <- glmmTMB(rate ~ length + I(length^2) + Tag_area + Release_year + julian_recapture_std, family=gaussian, data=Data_mackerel_use_Ireland_select)
      m15 <- glmmTMB(rate ~ I(log(length)) + I(log(length)^2) + Tag_area + Release_year+ I(julian_recapture_std), family=gaussian, data=Data_mackerel_use_Ireland_select)
      m16 <- glmmTMB(rate ~ I(log(length)) + I(log(length)^2) + Tag_area + Release_year+ julian_recapture_std, family=Gamma(link = "log"), data=Data_mackerel_use_Ireland_select)

      m15b <- glmmTMB(rate ~ length + Tag_area*Release_timing + julian_recapture_std:Release_year, family=Gamma(link = "log"), data=Data_mackerel_use_Ireland_select)

      m16 <- glmmTMB(log_rate ~ scale(length) + Tag_area*Release_timing + Release_year + scale(julian_recapture_std), family=gaussian, data=Data_mackerel_use_Ireland_select)

      AICtab(m1,m2,m3,m4,m5,m6,m7,m8,
             m0,m01,m02,m03,m04,m05,m10,m11,m12,m13)
      bestm <- m00
      plot(bestm, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_select <- simulateResiduals(fittedModel = bestm, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_select, rank=T, quantreg = TRUE)
      testResiduals(sim_select)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Ireland_select$Length, sim_select$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Ireland_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Ireland_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Ireland_select$Latitude, sim_select$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Ireland_select$Longitude, sim_select$scaledResiduals, main = "Longitude")

      ## Using duration std1
      # mmm0_Ireland <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm01_Ireland <- gam(log_rate_annual1 ~ s(length) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
      # mmm02_Ireland <- gam(log_rate_annual1 ~ s(length) + Category, family=gaussian, data=Data_mackerel_use_Ireland_select)
      # mmm03_Ireland <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_Ireland_select)
      # mmm04_Ireland <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
      # mmm05_Ireland <- gam(log_rate_annual1 ~ length + s(Longitude, Latitude, by= Release_year) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
      # mmm06_Ireland <- gam(log_rate_annual1 ~ length + s(Longitude, Latitude) + Release_year + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
      # AICtab(mmm0_Ireland, mmm01_Ireland,mmm02_Ireland,mmm03_Ireland,mmm04_Ireland,mmm05_Ireland,mmm06_Ireland)

      # mmm1_Ireland <- gam(log_rate_annual1 ~ s(length) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm2_Ireland <- gam(log_rate_annual1 ~ s(length) + Release_month + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm3_Ireland <- gam(log_rate_annual1 ~ s(length) + Release_month + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm4_Ireland <- gam(log_rate_annual1 ~ s(length, by= Release_year) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm5_Ireland <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm6_Ireland  <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland)
      # mmm7_Ireland  <- gam(log_rate_annual1 ~ s(length, by= Release_year) + s(Longitude, Latitude, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland)

      # AICtab(mmm0_Ireland, mmm1_Ireland,mmm2_Ireland,mmm3_Ireland,mmm4_Ireland,mmm5_Ireland,mmm6_Ireland,mmm7_Ireland)

      # best_m <- mmm05_Ireland
      # plot(best_m, all.terms=TRUE, residuals=TRUE, pages=1)
      # sim_std1_Ireland <- simulateResiduals(fittedModel = best_m, n = 1000, integerResponse = FALSE, plot=FALSE)
      # plot(sim_std1_Ireland, rank=T, quantreg = TRUE)
      # testResiduals(sim_std1_Ireland)
      # par(mfrow=c(3,2))
      # plotResiduals(Data_mackerel_use_Ireland$Length, sim_std1_Ireland$scaledResiduals, main = "length")
      # plotResiduals(Data_mackerel_use_Ireland$Release_month, sim_std1_Ireland$scaledResiduals, main = "Release_month")
      # plotResiduals(Data_mackerel_use_Ireland$Release_year, sim_std1_Ireland$scaledResiduals, main = "Release_year")
      # plotResiduals(Data_mackerel_use_Ireland$Latitude, sim_std1_Ireland$scaledResiduals, main = "Latitude")
      # plotResiduals(Data_mackerel_use_Ireland$Longitude, sim_std1_Ireland$scaledResiduals, main = "Longitude")

      ## Conclusion:
      ## Residual pattern does not look good


      with(Data_mackerel_use_Ireland_select, plot(length, log_rate))
      with(subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland"), plot(length, log_rate))
      with(subset(Data_mackerel_use_Ireland_select, Tag_area=="South_Ireland"), plot(length, log_rate))
      with(subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland"), plot(date, log_rate))
      with(subset(Data_mackerel_use_Ireland_select, Tag_area=="South_Ireland"), plot(date, log_rate))
      ggplot(Data_mackerel_use_Ireland_select, aes(x=date, y=log_rate, col=Release_year)) + geom_point() + facet_grid(~Tag_area)
      ggplot(Data_mackerel_use_Ireland_select, aes(x=date, y=log_rate, col=Release_year)) + geom_point() + facet_grid(~Tag_area) + geom_hline(yintercept = mean_diff_tag_area+ 8.7)+ geom_hline(yintercept =  8.7)

      brr <- subset(Data_mackerel_use_Ireland_select, Tag_area=="South_Ireland")
      brr$col <- ifelse(brr$log_rate<8.9,"red", ifelse(brr$log_rate>9.7, "blue", "purple"))
      ggmap(map_area) + geom_point(data=brr, aes(x=cLon, y=cLat, col=col), size=1)
      brr1 <- subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland")
      brr1$col <- ifelse(brr1$log_rate<(8.9-mean_diff_tag_area),"red", "purple")
      ggmap(map_area) + geom_point(data=brr1, aes(x=cLon, y=cLat, col=col), size=1)

      m0 <- gam(log_rate ~ s(length) + s(Longitude, Latitude), family=gaussian, data=subset(brr, col=="red"))
      simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)
      m0 <- gam(log_rate ~ s(length) + s(Longitude, Latitude), family=gaussian, data=subset(brr, col=="blue"))
      simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)
      m0 <- gam(log_rate ~ s(length) + s(Longitude, Latitude), family=gaussian, data=subset(brr, col=="purple"))
      simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)



  }


########## Preparing the final data from Ireland for the analysis

	Data_mackerel_use_Ireland1 <- Data_mackerel_use_Ireland
	Data_mackerel_use_Ireland1$Catch_month <- as.numeric(as.character(Data_mackerel_use_Ireland1$Catch_month))
	Data_mackerel_use_Ireland1$Catch_year <- as.numeric(as.character(Data_mackerel_use_Ireland1$Catch_year))
	Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(1,2)),'Catch_month'] <- Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(1,2)),'Catch_month'] + 12
	Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(13,14)),'Catch_year'] <- Data_mackerel_use_Ireland1[which(Data_mackerel_use_Ireland1$Catch_month %in% c(13,14)),'Catch_year'] - 1

	## Deriving data-frame for the 3 different migration cycles
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

	kfolding <- function(it, kfolds = 10, model_type="Iceland", inputdata=www, ...){

	  set.seed(it)
	  ERRORS <- c()
	  AICs <- c()
	  BICs <- c()
	  if(model_type == "Iceland") inputdata$y <- inputdata$toiceland
	  if(model_type == "Norway") inputdata$y <- inputdata$to_norway
	  if(model_type == "Northsea") inputdata$y <- inputdata$to_northsea
	  if(model_type == "Ireland") inputdata$y <- inputdata$to_ireland
	  www1 <- inputdata %>% ungroup() %>% fold(k = kfolds, method="n_rand")
	  dat_use <- www1

	  for (k in 1:kfolds){

	    testing <- www1[www1$.folds == k,]
	    training <- www1[-testing$ID,]

	    if (kfolds > 1) dat_use <- training

	    m0 <- (gam(y ~  1, family=binomial,
	               data = dat_use))
	    m01 <- (gam(y ~  -1 + Release_year, family=binomial,
	                data = dat_use))
	    m02 <- (gam(y ~  -1 + Release_year+ s(Latitude), family=binomial,
	                data = dat_use))
	    m03 <- (gam(y ~  -1 + Release_year+ s(Length), family=binomial,
	                data = dat_use))
	    m04 <- (gam(y ~  -1 + Release_year+ s(julian_recapture_std_scaled, k=3), family=binomial,
	                data = dat_use))
	    m05 <- (gam(y ~  -1 + Release_year+ s(Latitude) + s(julian_recapture_std_scaled, k=3), family=binomial,
	                data = dat_use))
	    m06 <- (gam(y ~  -1 + Release_year+ s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
	                data = dat_use))
	    m07 <- (gam(y ~  -1 + Release_year+ s(Latitude) + s(Length), family=binomial,
	                data = dat_use))
	    if (month <= 12) m1 <- (gam(y ~  -1 + s(Latitude) + Release_year  + s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
	               data = dat_use))
	    if (month > 12) m1 <- (gam(y ~  -1 + s(Latitude,k=3) + Release_year  + s(Length,k=3) + s(julian_recapture_std_scaled, k=3), family=binomial,
	               data = dat_use))
	    # m2 <- (gam(y ~  -1 + s(Latitude) + Release_year + s(Length, by=Release_year) + s(julian_recapture_std_scaled, k=3), family=binomial,
	    #            data = dat_use))
	    # m3 <- (gam(y ~  -1 + s(Latitude, by=Release_year) + Release_year + s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
	    #            data = dat_use))
	    # m4 <- (gam(y ~  -1 + s(Latitude) + Release_year  + s(Length) + s(julian_recapture_std_scaled, by=Release_year), family=binomial,
	    #            data = dat_use))
	    m2 <- (gam(y ~  -1 + s(Latitude) + Release_year + s(julian_recapture_std_scaled, by=length_bin), family=binomial,
	               data = dat_use))
	    m3 <- (gam(y ~  -1 + s(Latitude) + Release_year + s(Length, k=3) + s(julian_recapture_std_scaled, by=length_bin), family=binomial,
	               data = dat_use))
	    m4 <- (gam(y ~  -1 + s(Latitude, by=length_bin) + s(Length, k=3) + Release_year + s(julian_recapture_std_scaled, by=length_bin), family=binomial,
	               data = dat_use))
	    mse0 = mse01 = mse02 = mse03 = mse04 = mse05 = mse06 = mse07 = mse1 = mse2 = mse3 = mse4 = NA

	    pred0 <- predict(m0, type="response", newdata = testing)
	    pred01 <- predict(m01, type="response", newdata = testing)
	    pred02 <- predict(m02, type="response", newdata = testing)
	    pred03 <- predict(m03, type="response", newdata = testing)
	    pred04 <- predict(m04, type="response", newdata = testing)
	    pred05 <- predict(m05, type="response", newdata = testing)
	    pred06 <- predict(m06, type="response", newdata = testing)
	    pred07 <- predict(m07, type="response", newdata = testing)
	    pred1 <- predict(m1, type="response", newdata = testing)
	    pred2 <- predict(m2, type="response", newdata = testing)
	    pred3 <- predict(m3, type="response", newdata = testing)
	    pred4 <- predict(m4, type="response", newdata = testing)

	    obs <- tonum(testing$y)
	    mse0 <- mean(apply(pred0-obs, 1, function(x) sum(x^2)))
	    mse01 <- mean(apply(pred01-obs, 1, function(x) sum(x^2)))
	    mse02 <- mean(apply(pred02-obs, 1, function(x) sum(x^2)))
	    mse03 <- mean(apply(pred03-obs, 1, function(x) sum(x^2)))
	    mse04 <- mean(apply(pred04-obs, 1, function(x) sum(x^2)))
	    mse05 <- mean(apply(pred05-obs, 1, function(x) sum(x^2)))
	    mse06 <- mean(apply(pred06-obs, 1, function(x) sum(x^2)))
	    mse07 <- mean(apply(pred07-obs, 1, function(x) sum(x^2)))
	    mse1 <- mean(apply(pred1-obs, 1, function(x) sum(x^2)))
	    mse2 <- mean(apply(pred2-obs, 1, function(x) sum(x^2)))
	    mse3 <- mean(apply(pred3-obs, 1, function(x) sum(x^2)))
	    mse4 <- mean(apply(pred4-obs, 1, function(x) sum(x^2)))

	    nAICs <- AIC(m0,m01,m02,m03,m04,m05,m06,m07,m1,m2,m3,m4)
	    nBICs <- BIC(m0,m01,m02,m03,m04,m05,m06,m07,m1,m2,m3,m4)
	    error <- c(mse0,mse01,mse02,mse03,mse04,mse05,mse06,mse07,mse1, mse2, mse3, mse4)
	    ERRORS <- rbind(ERRORS, error)
	    AICs <- rbind(AICs, nAICs)
	    BICs <- rbind(BICs, nBICs)
	  }

	  RESULT <- list(m0=m0,m01=m01,m02=m02,m03=m03,m04=m04,m05=m05,m06=m06,m07=m07,m1=m1,m2=m2,m3=m3,m4=m4,ERRORS=ERRORS, AICs=AICs, BICs=BICs, data=inputdata)

	  if(kfolds == 1) return(RESULT)
	  if(kfolds > 1) return(ERRORS)
	}

	run_directionality <- function(month, data_origin="cycle1", model_selection = "none", models=NULL, sensitivity=FALSE){
  ## Doing the analysis of P(moving to ICES eco-regions)
  	  Limit_month <- month ## c(9,10,11)
      if (month == 9)   label <- "Sep30th"
      if (month == 10)  label <- "Oct31st"
      if (month == 11)  label <- "Nov30st"
      if (month == 12)  label <- "Dec31st"
      if (month == 13)  label <- "Jan31st"
      if (month == 14)  label <- "Feb28th"

      if (data_origin=="cycle1") { lag = "cycle 1"      ; year_lag = 0  }
      if (data_origin=="cycle2") { lag = "cycle 2"  ; year_lag = 1  }
      if (data_origin=="cycle3") { lag = "cycle 3"  ; year_lag = 2  }


  			# Full data
          if (data_origin=="cycle1") data = Data_mackerel_use_Ireland_select_origin
          if (data_origin=="cycle2") data = Data_mackerel_use_Ireland_select_origin_year1
          if (data_origin=="cycle3") data = Data_mackerel_use_Ireland_select_origin_year2
          data$month <- as.numeric(data$Catch_month)
  	      Data_mackerel_use_Ireland_select <- subset(data, Release_year %in% 2014:2020) %>% filter(month <= Limit_month)

  	    # Now calculating the intersection between recapture location and ICES ecoregions + tweaking because some points fall on land
  	      Data_mackerel_use_Ireland_select_sf <- Data_mackerel_use_Ireland_select %>% st_as_sf(coords = c("cLon","cLat"), crs=4326) %>%
  	        st_transform(new_proj)
  	      Data_mackerel_use_Ireland_select$X <- st_coordinates(Data_mackerel_use_Ireland_select_sf)[,1]
  	      Data_mackerel_use_Ireland_select$Y <- st_coordinates(Data_mackerel_use_Ireland_select_sf)[,2]
#
#   	      ggplot(ICESecoregion_cut_proj) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
#   	        geom_point(data=Data_mackerel_use_Ireland_select, aes(x=X, y=Y))
#   	      ggplot(ICESecoregion_cut) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
#   	        geom_point(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat))

  	      intersection <- st_intersects(Data_mackerel_use_Ireland_select_sf, ICESecoregion_proj, sparse=FALSE)
  	      area <- rep(NA, nrow(intersection))
  	      for (i in 1:nrow(intersection)) { if (length(which(intersection[i,]==TRUE)>0)) area[i]= which(intersection[i,]==TRUE) }
  	      # ggplot(ICESecoregion_cut) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
  	        # geom_point(data=Data_mackerel_use_Ireland_select[which(is.na(area)==TRUE),], aes(x=cLon, y=cLat))
  	      # a little adjustment of the area falling on land
  	      area[which(is.na(area == TRUE) & (Data_mackerel_use_Ireland_select[, 'cLon'] == -1.5))] = 9
  	      area[which(is.na(area == TRUE) & (Data_mackerel_use_Ireland_select[, 'cLon'] == 5.5))] = 11
          # points falling in faroes but re-allocated to norway and iceland
  	      area[which(area == 15 & Data_mackerel_use_Ireland_select$cLon > -5)] = 16
          area[which(area == 15 & Data_mackerel_use_Ireland_select$cLon < -5)] = 13
  	      Data_mackerel_use_Ireland_select$direction <- ICESecoregion_proj$Ecoregion[area]
  	      # ggplot(ICESecoregion_cut_proj) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
  	      #   geom_point(data=Data_mackerel_use_Ireland_select[which(is.na(area)==TRUE),], aes(x=X, y=Y))

  			Data_mackerel_use_Ireland_select$ID <- 1:nrow(Data_mackerel_use_Ireland_select)

  			# Rescaling parameters to ease interpretation
  			Data_mackerel_use_Ireland_select$julian_recapture_scaled <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std)
  			Data_mackerel_use_Ireland_select$julian_recapture_standardized <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std, center=FALSE)
  			Data_mackerel_use_Ireland_select$length_scaled <- scale(Data_mackerel_use_Ireland_select$Length)
  			Data_mackerel_use_Ireland_select$Latitude_scaled <- scale(Data_mackerel_use_Ireland_select$Latitude)
  			Data_mackerel_use_Ireland_select$Latitude2_scaled <- scale((Data_mackerel_use_Ireland_select$Latitude)^2)

  			Data_mackerel_use_Ireland_select <- Data_mackerel_use_Ireland_select %>%
  			    mutate(toiceland= ifelse(direction == "Icelandic Waters", 1, 0),
  			           to_norway = ifelse(direction == "Norwegian Sea", 1, 0),
  			           to_northsea = ifelse(direction == "Greater North Sea", 1, 0),
  			           to_ireland = ifelse(direction == "Celtic Seas", 1, 0)
        )

  			# Data_mackerel_use_Ireland_select %>% group_by(toiceland, Release_year) %>% summarize(n=n())
  			#
  			# Data_mackerel_use_Ireland_select %>% filter(toiceland == 1) %>% ggplot(aes(x=Latitude)) + geom_histogram()
  			#
  			#
  			# hist(Data_mackerel_use_Ireland_select$dist_toiceland, breaks=50, xlab="Distance from Norway coast (m)", main="")

  	## Analysis of the EW movement
  			www <- Data_mackerel_use_Ireland_select
  			www <- www %>% mutate(group = paste0(Release_year, "_", Release_timing),
  			                      group = as.factor(group),
  			                      Release_timing_fact = as.factor(Release_timing),
  			                      group = droplevels(group),
  			                      Release_year = droplevels(Release_year),
  			                      julian_recapture_std_scaled = scale(julian_recapture_std),
  			                      julian_recapture_scaled = scale(julian_recapture),
  			                      Latitude_scaled = scale(Latitude),
  			                      Length_scaled = scale(Length),
  			                      toiceland_bin = as.factor(toiceland),
  			                      tonorway_bin = as.factor(to_norway),
  			                      tonorthsea_bin = as.factor(to_northsea),
  			                      toireland_bin = as.factor(to_ireland)
  			                      )

  			### 10-fold cross-validation

  			library(groupdata2)
  			library(furrr)
  			library(progressr)


  			if (is.null(models) == TRUE) {
  			### Go to Iceland

          main_iceland = kfolding(it=1, kfolds=1, model_type="Iceland", inputdata=www)

          # CVs <- kfolding(it=1, kfolds = 10, model_type="Iceland")
          # apply(CVs, 2, mean, na.rm=T)
          # CVs <- kfolding(it=20, kfolds = 10, model_type="Iceland")
          # apply(CVs, 2, mean, na.rm=T)
#
#           main_iceland$AICs
#           concurvity(main_iceland$m1)
#           concurvity(main_iceland$m2)
#           concurvity(main_iceland$m3)
#           concurvity(main_iceland$m4)

          library(DHARMa)
          library(ggeffects)
          if (model_selection == "AIC") best_iceland <- main_iceland[[which.min(main_iceland$AICs$AIC[1:9])]]
          if (model_selection == "none") best_iceland <- main_iceland[[9]]
          # par(mfrow=c(2,2))
          # plot.gam(best_iceland, all.terms=T)
          simul <- simulateResiduals(best_iceland, plot=TRUE)
          par(mfrow=c(2,2))
          plot(www$Latitude, simul$scaledResiduals)
          plot(www$Length, simul$scaledResiduals)
          plot(www$julian_recapture_std, simul$scaledResiduals)
          plot(www$Release_year, simul$scaledResiduals)

          # mydf <- ggpredict(main_iceland$m1, terms = "Latitude", back.transform = FALSE)
          # ggplot(mydf, aes(x, predicted)) +
          #   geom_line() +
          #   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)

        ### Go to norway
          main_norway = kfolding(it=1, kfolds=1, model_type="Norway", inputdata=www)
  			  main_norway$AICs
  			  main_norway$BICs
  			  concurvity(main_norway$m1)
  			  concurvity(main_norway$m2)
  			  concurvity(main_norway$m3)
  			  concurvity(main_norway$m4)
  			  library(DHARMa)
  			  if (model_selection == "AIC") best_norway <- main_norway[[which.min(main_norway$AICs$AIC[1:9])]]
  			  if (model_selection == "none") best_norway <- main_norway[[9]]
  			  simul <- simulateResiduals(best_norway, plot=TRUE)
  			  # par(mfrow=c(2,2))
  			  # plot.gam(best_norway, all.terms=TRUE)
  			  par(mfrow=c(2,2))
  			  plot(www$Latitude, simul$scaledResiduals)
  			  plot(www$Length, simul$scaledResiduals)
  			  plot(www$julian_recapture_std, simul$scaledResiduals)
  			  plot(www$Release_year, simul$scaledResiduals)

  			  # CVs <- kfolding(it=1, kfolds = 10, model_type="Norway")
  			  # apply(CVs, 2, mean, na.rm=T)
  			  # CVs <- kfolding(it=20, kfolds = 10, model_type="Norway")
  			  # apply(CVs, 2, mean, na.rm=T)
  			  #


			  ### Go to north sea

  			  main_northsea = kfolding(it=1, kfolds=1, model_type="Northsea", inputdata=www)
  			  main_northsea$AICs
  			  main_northsea$BICs
  			  concurvity(main_northsea$m1)
  			  concurvity(main_northsea$m2)
  			  concurvity(main_northsea$m3)
  			  concurvity(main_northsea$m4)
  			  library(DHARMa)
  			  if (model_selection == "AIC") best_northsea <- main_northsea[[which.min(main_northsea$AICs$AIC[1:9])]]
  			  if (model_selection == "none") best_northsea <- main_northsea[[9]]
  			  simul <- simulateResiduals(best_northsea, plot=TRUE)
  			  # par(mfrow=c(2,2))
  			  # plot.gam(best_northsea, all.terms=TRUE)
  			  par(mfrow=c(2,2))
  			  plot(www$Latitude, simul$scaledResiduals)
  			  plot(www$Length, simul$scaledResiduals)
  			  plot(www$julian_recapture_std, simul$scaledResiduals)
  			  plot(www$Release_year, simul$scaledResiduals)

  			  # CVs <- kfolding(it=1, kfolds = 10, model_type="Northsea")
  			  # apply(CVs, 2, mean, na.rm=T)
  			  # CVs <- kfolding(it=20, kfolds = 10, model_type="Northsea")
  			  # apply(CVs, 2, mean, na.rm=T)

  			### Go to Ireland

  			  main_ireland = kfolding(it=1, kfolds=1, model_type="Ireland", inputdata=www)
  			  main_ireland$AICs
  			  # main_ireland$BIC
  			  # concurvity(main_ireland$m1)
  			  # concurvity(main_ireland$m2)
  			  # concurvity(main_ireland$m3)
  			  # concurvity(main_ireland$m4)
  			  library(DHARMa)
  			  if (model_selection == "AIC") best_ireland <- main_ireland[[which.min(main_ireland$AICs$AIC[1:9])]]
  			  if (model_selection == "none") best_ireland <- main_ireland[[9]]
  			  simul <- simulateResiduals(best_ireland, plot=TRUE)
  			  # par(mfrow=c(2,2))
  			  # plot.gam(best_ireland, all.terms=TRUE)
  			  par(mfrow=c(2,2))
  			  plot(www$Latitude, simul$scaledResiduals)
  			  plot(www$Length, simul$scaledResiduals)
  			  plot(www$julian_recapture_std, simul$scaledResiduals)
  			  plot(www$Release_year, simul$scaledResiduals)

  			  # CVs <- kfolding(it=1, kfolds = 10, model_type="Northsea")
  			  # apply(CVs, 2, mean, na.rm=T)
  			  # CVs <- kfolding(it=20, kfolds = 10, model_type="Northsea")
  			  # apply(CVs, 2, mean, na.rm=T)

  		}

  			if (is.null(models) == FALSE) {
  			  main_iceland = models$main_iceland
  			  main_ireland = models$main_ireland
  			  main_norway = models$main_norway
  			  main_northsea = models$main_northsea

  			  if (model_selection == "none"){
  			    best_iceland = main_iceland[[9]]
  			    best_ireland = main_ireland[[8]]
  			    best_norway = main_norway[[9]]
  			    best_northsea = main_northsea[[9]]
  			  }
  			  if (model_selection == "AIC"){
  			    best_iceland = main_iceland[[which.min(main_iceland$AICs$AIC[1:9])]]
  			    best_ireland = main_ireland[[which.min(main_ireland$AICs$AIC[1:9])]]
  			    best_norway = main_norway[[which.min(main_norway$AICs$AIC[1:9])]]
  			    best_northsea = main_northsea[[c(1,2,3,4,8)[which.min(main_ireland$AICs$AIC[c(1,2,3,4,8)])]]]
  			  }
  			}



      ### creating figures now
    		# some plotting configurations (either labels, or plot itself)
   			  int_breaks <- function(x, n = 5) {
    			  l <- pretty(x, n)
    			  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
    			}
  			  func1space <- function(x) paste0(" ", x)
  			  func2space <- function(x) paste0("  ", x)
  			  ggctr0 <- theme(axis.title = element_text(size=13),
  			                  axis.text = element_text(size=10),
  			                  axis.title.y=element_blank(),
  			                  strip.text.x = element_text(size = 10),
  			                  axis.text.y = element_text(angle =90, hjust=0.5))
  			  ggctr1 <- theme(axis.title = element_text(size=13),
  			                  axis.text = element_text(size=10),
  			                  strip.text.x = element_text(size = 10),
  			                  axis.text.y = element_text(angle =90, hjust=0.5))
  			  ggctr2 <- theme(axis.title = element_text(size=13),
  			                  axis.text = element_text(size=10),
  			                  axis.title.x=element_blank(),
  			                  strip.text.x = element_text(size = 10),
  			                  axis.text.y = element_text(angle =90, hjust=0.5))
  			  ggctr4 <- theme(axis.title = element_text(size=13),
  			                  axis.text = element_text(size=10),
  			                  axis.title.x=element_blank(),
  			                  axis.title.y=element_blank(),
  			                  strip.text.x = element_text(size = 10),
  			                  axis.text.y = element_text(angle =90, hjust=0.5))
  			  ggctrl3 <- theme(axis.title = element_text(size=13),
  			                   axis.text = element_text(size=10),
  			                   axis.text.x = element_text(size=10),
  			                   strip.text.x = element_text(size = 10),
  			                   axis.text.y = element_text(angle =90, hjust=0.5),
  			                   plot.tag.position = c(0.55, 0.07),
  			                   plot.tag = element_text(size=13))

  			 # determining the total sample size (with number of observation in the region in parenthesis)
  			  nsamp_IS <- www %>% group_by(toiceland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toiceland_bin == 1)
  			    nsamp_NO <- www %>% group_by(tonorway_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorway_bin == 1)
  			    nsamp_NS <- www %>% group_by(tonorthsea_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorthsea_bin == 1)
  			  if (Limit_month > 12 & sensitivity == FALSE){
  		      nsamp_IR <- www %>% group_by(toireland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toireland_bin == 1)
  			  }

  			# determinig the % deviance explained
  			    dev_IS <- data.frame(label=paste0("dev expl: ", round(summary(best_iceland)$dev.expl*100,0), "%"))
  			    dev_NO <- data.frame(label=paste0("dev expl: ", round(summary(best_norway)$dev.expl*100,0), "%"))
  			    dev_NS <- data.frame(label=paste0("dev expl: ", round(summary(best_northsea)$dev.expl*100,0), "%"))
  			    if (Limit_month > 12 & sensitivity == FALSE) dev_IR <- data.frame(label=paste0("dev expl: ", round(summary(best_ireland)$dev.expl*100,0), "%"))

        # Latitude effect
  			if (length(grep("Latitude", best_iceland$call))>0) {
    			pp_IS <- visreg::visreg(fit=best_iceland, xvar="Latitude", plot=FALSE, data=main_iceland$data)
    			p0_IS <- ggplot(pp_IS$fit, aes(x=Latitude, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="", x="Latitude (º)") +
    			  theme(axis.title = element_text(size=13),
    			        axis.text = element_text(size=10),
    			        axis.title.x=element_blank(),
    			        strip.text.x = element_text(size = 10),
    			        axis.text.y = element_text(angle =90, hjust=0.5))
    			pp0_IS <- p0_IS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==0), aes(x=Latitude)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==1), aes(x=Latitude), col="red") #+
    			  #geom_text(data=nsamp_IS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			} else { pp0_IS <- ggplot() + theme_void() }
        if (length(grep("Latitude", best_norway$call))>0) {
          pp_NO <- visreg::visreg(best_norway, "Latitude", plot=FALSE, data=main_norway$data)
    			p0_NO <- ggplot(pp_NO$fit, aes(x=Latitude, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="Marginal effects", x="Latitude (º)") +
    			  theme(axis.title = element_text(size=13),
    			        axis.text = element_text(size=10),
    			        axis.title.x=element_blank(),
    			        strip.text.x = element_text(size = 10),
    			        axis.text.y = element_text(angle =90, hjust=0.5))#+
    			 # scale_y_continuous(labels=func1space)
    			pp0_NO <- p0_NO + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==0), aes(x=Latitude)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==1), aes(x=Latitude), col="red") #+
    			  #geom_text(data=nsamp_NO, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
        } else { pp0_NO <- ggplot() + theme_void() }
  			if (length(grep("Latitude", best_northsea$call))>0) {
    			pp_NS <- visreg::visreg(best_northsea, "Latitude", plot=FALSE, data=main_northsea$data)
    			p0_NS <- ggplot(pp_NS$fit, aes(x=Latitude, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="", x="Latitude (º)")
    			if( Limit_month <= 12 | sensitivity == TRUE) p0_NS <- p0_NS + ggctr1
    			if( Limit_month > 12 & sensitivity == FALSE) p0_NS <- p0_NS + ggctr2
    			pp0_NS <- p0_NS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=Latitude)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=Latitude), col="red") #+
    			  #geom_text(data=nsamp_NS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			} else {
  			  pp0_NS <- ggplot() + theme_void()
  			  if ( Limit_month <= 12 | sensitivity == TRUE) pp0_NS <- pp0_NS + labs(tag="Latitude (º)") + ggctrl3
  			  if ( Limit_month > 12 & sensitivity == FALSE) pp0_NS <- pp0_NS + ggctrl4
  			}
  			 if ( Limit_month > 12 & sensitivity == FALSE){
  			   if (length(grep("Latitude", best_ireland$call))>0) {
  			     pp_IR <- visreg::visreg(best_ireland, "Latitude", plot=FALSE, data=main_ireland$data)
  			     p0_IR <- ggplot(pp_IR$fit, aes(x=Latitude, y=visregFit)) + geom_line() +
  			       geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
  			       theme_bw() +
  			       labs(y="", x="Latitude (º)")
  			     if (data_origin!="cycle3") p0_IR <- p0_IR + ggctr1 #+ scale_y_continuous(breaks=int_breaks2, labels=func1space)
  			     if (data_origin=="cycle3") p0_IR <- p0_IR + ggctr1 #+ scale_y_continuous(breaks=int_breaks2, labels=func2space)
  			     pp0_IR <- p0_IR + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toireland_bin==0), aes(x=Latitude)) +
  			       geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toireland_bin==1), aes(x=Latitude), col="red") #+
  			     #geom_text(data=nsamp_IR, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			   } else {
  			     pp0_IR <- ggplot() + theme_void() + labs(tag="Latitude (º)") +
  			       theme(axis.title = element_text(size=13),
  			             axis.text = element_text(size=10),
  			             axis.text.x = element_text(size=10),
  			             strip.text.x = element_text(size = 10),
  			             plot.tag.position = c(0.55, 0.07),
  			             plot.tag = element_text(size=13))
   			   }
  			 }

        # Length effect
  			if (length(grep("Length", best_iceland$call))>0) {
    			pp_IS <- visreg::visreg(best_iceland, "Length", plot=FALSE, data=main_iceland$data)
    			p0_IS <- ggplot(pp_IS$fit, aes(x=Length, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="Marginal effect", x="Length (cm)") +
    			  theme(axis.title = element_text(size=13),
    			        axis.text = element_text(size=10),
    			        axis.title.x=element_blank(),
    			        axis.title.y=element_blank(),
    			        axis.text.y = element_text(angle =90, hjust=0.5),
    			        strip.text.x = element_text(size = 10))
    			pp1_IS <- p0_IS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==0), aes(x=Length)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==1), aes(x=Length), col="red")
  			} else { pp1_IS <- ggplot() + theme_void() }
        if (length(grep("Length", best_norway$call))>0) {
    			pp_NO <- visreg::visreg(best_norway, "Length", plot=FALSE, data=main_norway$data)
    			p0_NO <- ggplot(pp_NO$fit, aes(x=Length, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="Marginal effect", x="Length (cm)") +
    			  theme(axis.title = element_text(size=13),
    			        axis.text = element_text(size=10),
    			        axis.title.y=element_blank(),
    			        axis.title.x=element_blank(),
    			        axis.text.y = element_text(angle =90, hjust=0.5),
    			        strip.text.x = element_text(size = 10)) #+
    			  #scale_y_continuous(breaks=int_breaks)
    			pp1_NO <- p0_NO + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==0), aes(x=Length)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==1), aes(x=Length), col="red")
        } else { pp1_NO <- ggplot() + theme_void() }
  			if (length(grep("Length", best_northsea$call))>0) {
    			pp_NS <- visreg::visreg(best_northsea, "Length", plot=FALSE, data=main_northsea$data)
    			p0_NS <- ggplot(pp_NS$fit, aes(x=Length, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="Marginal effect", x="Length (cm)")
    			if( Limit_month <= 12 | sensitivity == TRUE) p0_NS <- p0_NS + ggctr0 #+
    			  #scale_y_continuous(labels=func1space)
    			if( Limit_month > 12 & sensitivity == FALSE) p0_NS <- p0_NS + ggctr0 + theme(axis.title.x = element_blank())
    			pp1_NS <- p0_NS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=Length)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=Length), col="red")
  			} else {
  			  pp1_NS <- ggplot() + theme_void()
  			  if ( Limit_month <= 12 | sensitivity == TRUE) pp1_NS <- pp1_NS + labs(tag="Length (cm)") + ggctrl3
  			  if ( Limit_month > 12 & sensitivity == FALSE) pp1_NS <- pp1_NS + ggctrl4
         }
        if ( Limit_month > 12 & sensitivity == FALSE){
          if (length(grep("Length", best_ireland$call))>0) {
            pp_IR <- visreg::visreg(best_ireland, "Length", plot=FALSE, data=main_ireland$data)
            p1_IR <- ggplot(pp_IR$fit, aes(x=Length, y=visregFit)) + geom_line() +
              geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
              theme_bw() +
              labs(y="", x="Length (cm)")
            p1_IR <- p1_IR + ggctr0 # + scale_y_continuous(breaks=int_breaks, labels=func1space)
            pp1_IR <- p1_IR + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toireland_bin==0), aes(x=Length)) +
              geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toireland_bin==1), aes(x=Length), col="red")
          } else {
            pp1_IR <- ggplot() + theme_void() + labs(tag="Length (cm)") +
              theme(axis.title = element_text(size=13),
                    axis.text = element_text(size=10),
                    axis.text.x = element_text(size=10),
                    strip.text.x = element_text(size = 10),
                    plot.tag.position = c(0.55, 0.07),
                    plot.tag = element_text(size=13))
          }
        }

        # Recapture date effect
  			if (length(grep("recapture", best_iceland$call))>0) {
    			pp_IS <- visreg::visreg(best_iceland, "julian_recapture_std_scaled", plot=FALSE, data=main_iceland$data)
    			p0_IS <- ggplot(pp_IS$fit, aes(x=julian_recapture_std_scaled, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="Marginal effect", x="Recapture date std") +
    			  theme(axis.title = element_text(size=13),
    			        axis.text = element_text(size=10),
    			        axis.title.y=element_blank(),
    			        axis.title.x=element_blank(),
    			        axis.text.y = element_text(angle =90, hjust=0.5),
    			        strip.text.x = element_text(size = 10))
    			pp2_IS <- p0_IS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==0), aes(x=julian_recapture_std_scaled)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==1), aes(x=julian_recapture_std_scaled), col="red") #+
    			  #geom_text(data=nsamp_IS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			} else { pp2_IS <- ggplot() + theme_void() }
        if (length(grep("recapture", best_norway$call))>0) {
    			pp_NO <- visreg::visreg(best_norway, "julian_recapture_std_scaled", plot=FALSE, data=main_norway$data)
    			p0_NO <- ggplot(pp_NO$fit, aes(x=julian_recapture_std_scaled, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw() +
    			  labs(y="Marginal effect", x="Recapture date std") +
    			  theme(axis.title = element_text(size=13),
    			        axis.text = element_text(size=10),
    			        axis.title.y=element_blank(),
    			        axis.title.x=element_blank(),
    			        axis.text.y = element_text(angle =90, hjust=0.5),
    			        strip.text.x = element_text(size = 10))
    			pp2_NO <- p0_NO + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==0), aes(x=julian_recapture_std_scaled)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==1), aes(x=julian_recapture_std_scaled), col="red")# +
    			  #geom_text(data=nsamp_NO, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
        } else { pp2_NO <- ggplot() + theme_void() }
  			if (length(grep("recapture", best_northsea$call))>0) {
    			pp_NS <- visreg::visreg(best_northsea, "julian_recapture_std_scaled", plot=FALSE, data=main_northsea$data)
    			p0_NS <- ggplot(pp_NS$fit, aes(x=julian_recapture_std_scaled, y=visregFit)) + geom_line() +
    			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
    			   theme_bw()+ labs(y="Marginal effects", x="Recapture date std")
    			if( Limit_month <= 12 | sensitivity == TRUE) p0_NS <- p0_NS + ggctr0 #+
    			 # scale_y_continuous(labels=func1space)
    			if( Limit_month > 12 & sensitivity == FALSE) p0_NS <- p0_NS + ggctr0 + theme(axis.title.x = element_blank())
    			pp2_NS <- p0_NS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=julian_recapture_std_scaled)) +
    			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=julian_recapture_std_scaled), col="red")# +
    			  #geom_text(data=nsamp_NS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			} else {
  			  pp2_NS <- ggplot() + theme_void()
  			  if ( Limit_month <= 12 | sensitivity == TRUE) pp2_NS <- pp2_NS + labs(tag="Recapture date std") + ggctrl3
  			  if ( Limit_month > 12 & sensitivity == FALSE) pp2_NS <- pp2_NS + ggctrl4
  			}
        if ( Limit_month > 12 & sensitivity == FALSE){
          if (length(grep("julian_recapture_std_scaled", best_ireland$call))>0) {
            pp_IR <- visreg::visreg(best_ireland, "julian_recapture_std_scaled", plot=FALSE, data=main_ireland$data)
            p1_IR <- ggplot(pp_IR$fit, aes(x=julian_recapture_std_scaled, y=visregFit)) + geom_line() +
              geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
              theme_bw() +
              labs(y="", x="Recapture date std")
            p1_IR <- p1_IR + ggctr0
            pp2_IR <- p1_IR + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toireland_bin==0), aes(x=julian_recapture_std_scaled)) +
              geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toireland_bin==1), aes(x=julian_recapture_std_scaled), col="red") #+
            #geom_text(data=nsamp_IR, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
          } else {
            pp2_IR <- ggplot() + theme_void() + labs(tag="Recapture date std") +
              theme(axis.title = element_text(size=13),
                    axis.text = element_text(size=10),
                    axis.text.x = element_text(size=10),
                    strip.text.x = element_text(size = 10),
                    plot.tag.position = c(0.55, 0.07),
                    plot.tag = element_text(size=13))
          }
        }

        # year effect
  			  nsamp_year_IS <- www %>% group_by(Release_year, toiceland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Release_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toiceland_bin == 1)
   			  nsamp_year_NO <- www %>% group_by(Release_year, tonorway_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Release_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorway_bin == 1)
   			  nsamp_year_NS <- www %>% group_by(Release_year, tonorthsea_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Release_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorthsea_bin == 1)
   			  nsamp_year_IR <- www %>% group_by(Release_year, toireland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Release_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toireland_bin == 1)

   			  if (length(grep("Release_year", best_iceland$call))>0) {
   			    pp <- visreg::visreg(best_iceland, "Release_year", plot=FALSE, data=main_iceland$data)
      			pp$fit <- pp$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
   			    pp3_IS <- ggplot(pp$fit, aes(x=Release_year_fct, y=visregFit)) + geom_point() +
      			  geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col=grey(0.5), width=0.4) +
      			  theme_bw() +scale_x_discrete(drop=FALSE) +
      			  ylab("") +
      			  xlab("Year") + labs(tag="P(Icelandic waters)") +
      			  theme(axis.title = element_text(size=13),
      			        axis.text = element_text(size=10),
      			        axis.title.y=element_blank(),
      			        axis.title.x=element_blank(),
      			        axis.text.x = element_text(size = 7),
      			        axis.text.y = element_text(angle =90, hjust=0.5),
      			        plot.margin = margin(0.2,1.2,0.25,0.2, "cm"),
      			        plot.tag.position = c(1.1, 0.5),
      			        plot.tag = element_text(angle=270, size=13)) +
      			  geom_text(data=nsamp_IS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3) #+
      			  #geom_text(data=dev_IS, aes(x=Inf, y=Inf, label=label), vjust=1.2, hjust=1.1, size=3)
   			  } else { pp3_IS <- ggplot() + theme_void() }
   			  if (length(grep("Release_year", best_norway$call))>0) {
   			    pp <- visreg::visreg(best_norway, "Release_year", plot=FALSE, data=main_norway$data)
   			    pp$fit <- pp$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
   			    pp3_NO <- ggplot(pp$fit, aes(x=Release_year_fct, y=visregFit)) + geom_point() +
      			  geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col=grey(0.5), width=0.4) +
      			  theme_bw() + scale_x_discrete(drop=FALSE) +
      			  ylab("") +
      			  xlab("Year") + labs(tag="P(Norwegian Sea)") +
      			  theme(axis.title = element_text(size=13),
      			        axis.text = element_text(size=10),
      			        axis.title.y=element_blank(),
      			        axis.title.x=element_blank(),
      			        axis.text.x = element_text(size = 7),
      			        axis.text.y = element_text(angle =90, hjust=0.5),
      			        plot.margin = margin(0.2,1.2,0.25,0.2, "cm"),
      			        plot.tag.position = c(1.1, 0.5),
      			        plot.tag = element_text(angle=270, size=13))+
      			  #scale_y_continuous(labels=func1space) +
      			  geom_text(data=nsamp_NO, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)# +
      			  #geom_text(data=dev_NO, aes(x=Inf, y=Inf, label=label), vjust=1.2, hjust=1.1, size=3)
   			  } else { pp3_NO <- ggplot() + theme_void() }
   			  if (length(grep("Release_year", best_northsea$call))>0) {
   			    pp <- visreg::visreg(best_northsea, "Release_year", plot=FALSE, data=main_northsea$data)
   			    pp$fit <- pp$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
   			    pp3_NS <- ggplot(pp$fit, aes(x=Release_year_fct, y=visregFit)) + geom_point() +
       			  geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col=grey(0.5), width=0.4) +
      			  theme_bw() +scale_x_discrete(drop=FALSE) +
      			  ylab("") +
      			  xlab("Year") + labs(tag="P(Greater Northsea)") +
      			  theme(axis.title = element_text(size=13),
      			        axis.text = element_text(size=10),
      			        axis.title.y=element_blank(),
      			        axis.title.x=element_blank(),
      			        axis.text.x = element_text(size = 7),
      			        axis.text.y = element_text(angle =90, hjust=0.5),
      			        plot.margin = margin(0.2,1.2,0.25,0.2, "cm"),
      			        plot.tag.position = c(1.1, 0.5),
      			        plot.tag = element_text(angle=270, size=13))+
      			  #scale_y_continuous(breaks=int_breaks, labels=func2space)  +
      			  geom_text(data=nsamp_NS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3) #+
      			  #geom_text(data=dev_NS, aes(x=Inf, y=Inf, label=label), vjust=1.2, hjust=1.1, size=3)
   			  } else { pp3_NS <- ggplot() + theme_void() }
   			  if ( Limit_month > 12 & sensitivity == FALSE){
   			    if (length(grep("Release_year", best_ireland$call))>0) {
   			      pp_IR <- visreg::visreg(best_ireland, "Release_year", plot=FALSE, data=main_ireland$data)
   			      pp$fit <- pp$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
   			      pp3_IR <- ggplot(pp$fit, aes(x=Release_year_fct, y=visregFit)) + geom_point() +
   			        geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col=grey(0.5), width=0.4) +
   			        theme_bw() +scale_x_discrete(drop=FALSE) +
   			        ylab("") +
   			        xlab("") + labs(tag="P(Celtic Seas)") +
   			        theme(axis.title = element_text(size=13),
   			              axis.text = element_text(size=10),
   			              axis.title.y=element_blank(),
   			              axis.text.x = element_text(size = 7),
   			              axis.text.y = element_text(angle =90, hjust=0.5),
   			              plot.margin = margin(0.2,1.2,0.25,0.2, "cm"),
   			              plot.tag.position = c(1.1, 0.6),
   			              plot.tag = element_text(angle=270, size=13))+
   			        #scale_y_continuous(breaks=int_breaks, labels=func1space)+
   			        geom_text(data=nsamp_IR, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)# +
   			        #geom_text(data=dev_IR, aes(x=Inf, y=Inf, label=label), vjust=1.2, hjust=1.1, size=3)
   			    } else {
   			      pp3_IR <- ggplot() + theme_void() + labs(tag="Year") +
   			        theme(axis.title = element_text(size=13),
   			              axis.text = element_text(size=10),
   			              axis.text.x = element_text(size=10),
   			              strip.text.x = element_text(size = 10),
   			              plot.tag.position = c(0.55, 0.07),
   			              plot.tag = element_text(size=13))
   			    }
   			  }

   		### Now group all figures
        if (Limit_month <= 12 | sensitivity == TRUE){
     			ppp <- grid.arrange(pp0_IS,pp1_IS,pp2_IS,pp3_IS,
                              pp0_NO,pp1_NO,pp2_NO,pp3_NO,
                              pp0_NS,pp1_NS,pp2_NS,pp3_NS,
                              layout_matrix=matrix(c(1:12), byrow=T, ncol=4, nrow=3), widths=c(1.07,1,1,1.1), heights=c(1,1,1.07))
          ggsave(ppp, file=paste0(getwd(), "/MS/figs/Marginal_cutoffmonth", month, "_lag", data_origin, model_selection, ".pdf"),width=26, height=20, units="cm", dpi = 450)
        }
        if ( Limit_month > 12 & sensitivity == FALSE){
     			ppp <- grid.arrange(pp0_IS,pp1_IS,pp2_IS,pp3_IS,
                              pp0_NO,pp1_NO,pp2_NO,pp3_NO,
                              pp0_NS,pp1_NS,pp2_NS,pp3_NS,
                              pp0_IR,pp1_IR,pp2_IR,pp3_IR,
                              layout_matrix=matrix(c(1:16), byrow=T, ncol=4, nrow=4), widths=c(1.07,1,1,1.1), heights=c(1,1,1,1.07))
          ggsave(ppp, file=paste0(getwd(), "/MS/figs/Marginal_cutoffmonth", month, "_lag", data_origin, model_selection, ".pdf"),width=26, height=20, units="cm", dpi = 450)
     # 			ppp <- grid.arrange(pp0_IR,pp1_IR,pp3_IR,
     #                          layout_matrix=matrix(c(1:3), byrow=T, ncol=3, nrow=1), widths=c(1.07,1,1.1), heights=c(1))
     #      ggsave(ppp, file=paste0(getwd(), "/MS/figs/Marginal_cutoffmonth", month, "_lag", data_origin, model_selection, "IR.pdf"),width=26, height=7, units="cm", dpi = 450)
        }


        if (Limit_month <= 12 | sensitivity == TRUE) best_models <- list(main_iceland, main_norway, main_northsea)
        if ( Limit_month > 12 & sensitivity == FALSE) best_models <- list(main_iceland, main_norway, main_northsea, main_ireland, ppp)

      return(best_models)
  	}

	Resid_plot <- function(mod, dat, nsim=5000){
	  simul <- simulateResiduals(mod, plot=FALSE, n=nsim)
    par(mfrow=c(2,3), mar=c(4,3,1,0), oma=c(0,0,1,0))
    plotQQunif(simul)
    dat$residuals <- simul$scaledResiduals
    dat$Release_year <- droplevels(dat$Release_year)
    plot(dat$Length, dat$residuals, xlab="Fish bodysize (cm)", ylab="Scaled residuals");# lines(dat$Length, predict(gam(residuals ~ s(Length), data=dat)))
	  plot(dat$julian_recapture_std, dat$residuals, xlab="Recapture date", ylab="Scaled residuals"); #lines(dat$julian_recapture_std, predict(gam(residuals ~ s(julian_recapture_std), data=dat)))
	  plot(dat$Release_year, dat$residuals, xlab="Release year", ylab="Scaled residuals");
	  plot(dat$julian_release_std, dat$residuals, xlab="Release date", ylab="Scaled residuals"); #lines(dat$julian_release_std, predict(gam(residuals ~ s(julian_release_std), data=dat)))
	  # dat$resid <- simul$scaledResiduals
	  # dat$x <- simul$scaledResiduals
	  # d <- dat[order(dat$resid),]
	  # d$qq <- (qunif(ppoints(nrow(d))))
	  # p1 <- ggplot(d, aes(x= qq, y=resid)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw()
	  # p1 <- ggplot()
	}

	# Main run
  	Dec_0lag <- run_directionality(month=14, data_origin="cycle1", model_selection = "none", models=NULL)
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_0lag_IS.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_0lag[[1]][[9]], dat=Dec_0lag[[1]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_0lag_NO.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_0lag[[2]][[9]], dat=Dec_0lag[[2]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_0lag_NS.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_0lag[[3]][[9]], dat=Dec_0lag[[3]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_0lag_IR.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_0lag[[4]][[9]], dat=Dec_0lag[[4]]$data)
  	dev.off()

  	Dec_1lag <- run_directionality(month=14, data_origin="cycle2", model_selection = "none", models=NULL)
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_1lag_IS.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_1lag[[1]][[9]], dat=Dec_1lag[[1]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_1lag_NO.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_1lag[[2]][[9]], dat=Dec_1lag[[2]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_1lag_NS.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_1lag[[3]][[9]], dat=Dec_1lag[[3]]$data, nsim=5000)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_1lag_IR.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_1lag[[4]][[9]], dat=Dec_1lag[[4]]$data, nsim=5000)
  	dev.off()

  	Dec_2lag <- run_directionality(month=14, data_origin="cycle3", model_selection = "none", models=NULL)
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_2lag_IS.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_2lag[[1]][[9]], dat=Dec_2lag[[1]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_2lag_NO.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_2lag[[2]][[9]], dat=Dec_2lag[[2]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_2lag_NS.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_2lag[[3]][[9]], dat=Dec_2lag[[3]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/SuppMatDec_2lag_IR.png"), res=400, width=12, height=10, units="cm")
  	Resid_plot(mod=Dec_2lag[[4]][[9]], dat=Dec_2lag[[4]]$data)
  	dev.off()

  # main models based on model selection
  	Dec_1lag_AIC <- run_directionality(month=14, data_origin="cycle2", model_selection = "AIC", models=NULL)
  	Dec_2lag_AIC <- run_directionality(month=14, data_origin="cycle3", model_selection = "AIC", models=NULL)
  	Dec_0lag_AIC <- run_directionality(month=14, data_origin="cycle1", model_selection = "AIC", models=NULL)


	# Sensitivity
  	# Nov_0lag <- run_directionality(month=11, data_origin="cycle1", model_selection = "none", models=NULL)
  	# Nov_1lag <- run_directionality(month=11, data_origin="cycle2", model_selection = "none", models=NULL)
  	# Nov_2lag <- run_directionality(month=11, data_origin="cycle3", model_selection = "none", models=NULL)
  	Jan_0lagb <- run_directionality(month=13, data_origin="cycle1", model_selection = "none", models=NULL, sensitivity = FALSE)
  	Jan_1lagb <- run_directionality(month=13, data_origin="cycle2", model_selection = "none", models=NULL, sensitivity = FALSE)
  	Jan_2lagb <- run_directionality(month=13, data_origin="cycle3", model_selection = "none", models=NULL, sensitivity = FALSE)


##### Now doing some basic boxplots, etc for showing variability and tendency on mackerel length and latitude of release vs other covariates

  	create_data <- function(month=14, data_origin="cycle1"){
  	  Limit_month <- month ## c(9,10,11)
  	  if (month == 9)   label <- "Sep30th"
  	  if (month == 10)  label <- "Oct31st"
  	  if (month == 11)  label <- "Nov30st"
  	  if (month == 12)  label <- "Dec31st"
  	  if (month == 13)  label <- "Jan31st"
  	  if (month == 14)  label <- "Feb28th"

  	  if (data_origin=="cycle1") { lag = "cycle 1"      ; year_lag = 0  }
  	  if (data_origin=="cycle2") { lag = "cycle 2"  ; year_lag = 1  }
  	  if (data_origin=="cycle3") { lag = "cycle 3"  ; year_lag = 2  }


  	  # Full data
  	  if (data_origin=="cycle1") data = Data_mackerel_use_Ireland_select_origin
  	  if (data_origin=="cycle2") data = Data_mackerel_use_Ireland_select_origin_year1
  	  if (data_origin=="cycle3") data = Data_mackerel_use_Ireland_select_origin_year2
  	  data$month <- as.numeric(data$Catch_month)
  	  Data_mackerel_use_Ireland_select <- subset(data, Release_year %in% 2014:2020) %>% filter(month <= Limit_month)

  	  # Now calculating the intersection between recapture location and ICES ecoregions + tweaking because some points fall on land
  	  Data_mackerel_use_Ireland_select_sf <- Data_mackerel_use_Ireland_select %>% st_as_sf(coords = c("cLon","cLat"), crs=4326) %>%
  	    st_transform(new_proj)
  	  Data_mackerel_use_Ireland_select$X <- st_coordinates(Data_mackerel_use_Ireland_select_sf)[,1]
  	  Data_mackerel_use_Ireland_select$Y <- st_coordinates(Data_mackerel_use_Ireland_select_sf)[,2]
  	  #
  	  #   	      ggplot(ICESecoregion_cut_proj) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
  	  #   	        geom_point(data=Data_mackerel_use_Ireland_select, aes(x=X, y=Y))
  	  #   	      ggplot(ICESecoregion_cut) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
  	  #   	        geom_point(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat))

  	  intersection <- st_intersects(Data_mackerel_use_Ireland_select_sf, ICESecoregion_proj, sparse=FALSE)
  	  area <- rep(NA, nrow(intersection))
  	  for (i in 1:nrow(intersection)) { if (length(which(intersection[i,]==TRUE)>0)) area[i]= which(intersection[i,]==TRUE) }
  	  # ggplot(ICESecoregion_cut) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
  	  # geom_point(data=Data_mackerel_use_Ireland_select[which(is.na(area)==TRUE),], aes(x=cLon, y=cLat))
  	  # a little adjustment of the area falling on land
  	  area[which(is.na(area == TRUE) & (Data_mackerel_use_Ireland_select[, 'cLon'] == -1.5))] = 9
  	  area[which(is.na(area == TRUE) & (Data_mackerel_use_Ireland_select[, 'cLon'] == 5.5))] = 11
  	  # points falling in faroes but re-allocated to norway and iceland
  	  area[which(area == 15 & Data_mackerel_use_Ireland_select$cLon > -5)] = 16
  	  area[which(area == 15 & Data_mackerel_use_Ireland_select$cLon < -5)] = 13
  	  Data_mackerel_use_Ireland_select$direction <- ICESecoregion_proj$Ecoregion[area]
  	  # ggplot(ICESecoregion_cut_proj) + geom_sf(aes(geometry = geometry, fill=Ecoregion)) +
  	  #   geom_point(data=Data_mackerel_use_Ireland_select[which(is.na(area)==TRUE),], aes(x=X, y=Y))

  	  Data_mackerel_use_Ireland_select$ID <- 1:nrow(Data_mackerel_use_Ireland_select)

  	  # Rescaling parameters to ease interpretation
  	  Data_mackerel_use_Ireland_select$julian_recapture_scaled <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std)
  	  Data_mackerel_use_Ireland_select$julian_recapture_standardized <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std, center=FALSE)
  	  Data_mackerel_use_Ireland_select$length_scaled <- scale(Data_mackerel_use_Ireland_select$Length)
  	  Data_mackerel_use_Ireland_select$Latitude_scaled <- scale(Data_mackerel_use_Ireland_select$Latitude)
  	  Data_mackerel_use_Ireland_select$Latitude2_scaled <- scale((Data_mackerel_use_Ireland_select$Latitude)^2)

  	  Data_mackerel_use_Ireland_select <- Data_mackerel_use_Ireland_select %>%
  	    mutate(toiceland= ifelse(direction == "Icelandic Waters", 1, 0),
  	           to_norway = ifelse(direction == "Norwegian Sea", 1, 0),
  	           to_northsea = ifelse(direction == "Greater North Sea", 1, 0),
  	           to_ireland = ifelse(direction == "Celtic Seas", 1, 0)
  	    )

  	  # Data_mackerel_use_Ireland_select %>% group_by(toiceland, Release_year) %>% summarize(n=n())
  	  #
  	  # Data_mackerel_use_Ireland_select %>% filter(toiceland == 1) %>% ggplot(aes(x=Latitude)) + geom_histogram()
  	  #
  	  #
  	  # hist(Data_mackerel_use_Ireland_select$dist_toiceland, breaks=50, xlab="Distance from Norway coast (m)", main="")

  	  ## Analysis of the EW movement
  	  www <- Data_mackerel_use_Ireland_select
  	  www <- www %>% mutate(group = paste0(Release_year, "_", Release_timing),
  	                        group = as.factor(group),
  	                        Release_timing_fact = as.factor(Release_timing),
  	                        group = droplevels(group),
  	                        Release_year = droplevels(Release_year),
  	                        julian_recapture_std_scaled = scale(julian_recapture_std),
  	                        julian_recapture_scaled = scale(julian_recapture),
  	                        Latitude_scaled = scale(Latitude),
  	                        Length_scaled = scale(Length),
  	                        toiceland_bin = as.factor(toiceland),
  	                        tonorway_bin = as.factor(to_norway),
  	                        tonorthsea_bin = as.factor(to_northsea),
  	                        toireland_bin = as.factor(to_ireland)
  	  )  %>% ungroup() %>% dplyr::select(Release_year, Length, direction, Latitude, julian_release_std
    	                                  )
      www$cycle = lag

      return(www)
    }

    lag0data <- create_data(month=14, data_origin="cycle1")
    lag1data <- create_data(month=14, data_origin="cycle2")
    lag2data <- create_data(month=14, data_origin="cycle3")

    All_lags <- rbind(lag0data, lag1data, lag2data)
    All_lags <- All_lags %>% filter(direction != "Greenland Sea")

    library(corrplot)
    corrplot(cor(All_lags), method="number")

    glimpse(All_lags)
    library(ranger)
    All_lags <- All_lags
    All_lags_new <- All_lags %>% mutate(direction = "All")
    All_lags_full <- rbind(All_lags,All_lags_new ) %>% mutate(direction_short = as.factor(direction)) %>%
      mutate(direction_short = recode_factor(direction_short, "All"="All", "Icelandic Waters" = "ICW", "Norwegian Sea"="NWS", "Greater North Sea"="GNS", "Celtic Seas"="CES"),
             cycle = as.factor(cycle),
             Release_year = as.factor(Release_year))
    nsamp <- All_lags_full %>% group_by(direction_short,Release_year,cycle) %>% summarize(n = n()) %>%
      mutate(nsamp = paste0("n=", n), Latitude=50)
      # mutate(nsamp = paste0("n=", n), Latitude=ifelse(Release_year == 2020, 54, 60))
      # mutate(nsamp = ifelse(n<5, "*", ""), Latitude= 54)

    p1 <- All_lags_full %>% group_by(direction_short,Release_year,cycle) %>% mutate(mediany = median(Latitude)) %>%
      ggplot(., aes(y= Latitude, x = direction_short)) +
      geom_boxplot(aes(y= Latitude, x = direction_short), outlier.color=grey(0.8), outlier.shape = 2) +
      # geom_violin(scale ="area") + #, outlier.color=grey(0.8), outlier.shape = 2) +
      # geom_point(aes(x= direction_short, y=mediany), size=2) +
      geom_text(data=nsamp, aes(y= Latitude, x = direction_short, label=n), angle=0, size=2.5) +
      facet_grid(cycle ~ Release_year, scales="free_y") + theme_bw(18)+
      theme(axis.text.x = element_text(angle = 90, size=11)) + scale_y_continuous(limits = c(50, NA)) +
      ylab("Latitude (º)")
    ggsave(p1, file="MS/Figs/Latitude_distribution_boxplot.pdf", width=12, height=8, dpi=450)

    nsamp <- All_lags_full %>% group_by(direction_short,Release_year,cycle) %>% summarize(n = n()) %>%
      mutate(nsamp = paste0("n=", n), Length = 25)
      # mutate(nsamp = ifelse(n<5, "*", ""), Length= 42)
    p2 <- All_lags_full %>% group_by(direction_short,Release_year,cycle) %>% mutate(lengthy = median(Length)) %>%
      ggplot(.) +
      geom_boxplot(aes(y= Length, x = direction_short), outlier.color=grey(0.8), outlier.shape = 2) +
      # geom_violin(scale ="area", aes(y= Length, x = direction_short)) + #, outlier.color=grey(0.8), outlier.shape = 2) +
      # geom_point(aes(x= direction_short, y=lengthy), size=2) +
      geom_text(data=nsamp, aes(y= Length, x = direction_short, label=n), angle=0, size=2.5) +
      facet_grid(cycle ~ Release_year, scales="free_y") + theme_bw(18)+
      theme(axis.text.x = element_text(angle = 90, size=11)) + ylab("Length (cm)")
    ggsave(p2, file="MS/Figs/Length_distribution_boxplot.pdf", width=12, height=8, dpi=450)


    All_lags %>% group_by(direction_short,Release_year,cycle) %>%
      summarise(mean = mean(Length), ymin=quantile(Length, 0.025), ymax=quantile(Length, 0.975)) %>%
      ggplot(.) + geom_errorbar(aes(x=direction_short, ymin=ymin, ymax=ymax)) +
      geom_point(aes(x= direction_short, y=mean), size=2) +
      facet_grid(cycle ~ Release_year) + theme_bw(18)+ ylab("Length (cm)") +
      theme(axis.text.x = element_text(angle = 90, size=9))

    All_lags %>% group_by(direction_short,Release_year,cycle) %>%
      summarise(n = n(), mean = mean(Length), sd=sd(Length)) %>% mutate(se = sd/sqrt(n), ymin = mean - 1.96*se, ymax= mean + 1.96*se) %>%
      ggplot(.) + geom_errorbar(aes(x=direction_short, ymin=ymin, ymax=ymax)) +
      geom_point(aes(x= direction_short, y=mean), size=2) +
      facet_grid(cycle ~ Release_year) + theme_bw(18)+ ylab("Length (cm)") +
      theme(axis.text.x = element_text(angle = 90, size=9))


    glm1 <- ranger(Length ~ ., data=All_lags %>% dplyr::select(! c(Latitude, julian_release_std)), importance = "impurity")
  	glm1 <- ranger(Length ~ ., data=All_lags %>% dplyr::select(! c(Latitude)), importance = "impurity", mtry=3, num.trees=10000)
  	glm1 <- ranger(Length ~ ., data=All_lags %>% dplyr::select(! c(Latitude)), importance = "impurity", mtry=1, num.trees=10000)
  	glm1 <- ranger(Length ~ ., data=All_lags, importance = "impurity")
  	glm1b <- (lm(Length ~ cycle + direction * Release_year * julian_release_std, data=All_lags))
  	DHARMa::simulateResiduals(glm1b, plot=TRUE)
  	summary(glm1b)
  	anova(glm1b)
  	partial(glm1b, pred.var = "direction", plot=TRUE, rug=TRUE)
  	partial(glm1b, pred.var = "cycle", plot=TRUE, rug=TRUE)
  	partial(glm1b, pred.var = "Release_year", plot=TRUE, rug=TRUE)
  	partial(glm1b, pred.var = c("direction", "Release_year"), plot=TRUE, chull = TRUE)

  	library(pdp)
  	library(vip)
    vip(glm1, geom = "point", horizontal = TRUE)

    p1 <- partial(glm1, pred.var = "direction", plot=TRUE, rug=TRUE)
    p2 <- partial(glm1, pred.var = "cycle", plot=TRUE, rug=TRUE)  # Figure 2
    p3 <- partial(glm1, pred.var = "Release_year", plot=TRUE, rug=TRUE)  # Figure 2
    p4 <- partial(glm1, pred.var = "julian_release_std", plot=TRUE, rug=TRUE, smooth=TRUE)  # Figure 2
    # p5 <- partial(glm1, pred.var = "julian_recapture_std", plot=TRUE, rug=TRUE, smooth=TRUE)  # Figure 2
    # p6 <- partial(glm1, pred.var = c("julian_release_std", "julian_recapture_std"), plot=TRUE, chull = TRUE)
    # p4 <- glm1 %>% partial(pred.var = "julian_release_std") %>% autoplot(smooth = TRUE, ylab = "yhat") +
    #   theme_light()
    pp1 <- grid.arrange(p3,p1,p2, ncol=2)
    pp1 <- grid.arrange(p4, p3,p1,p2, ncol=2)
    pp1 <- grid.arrange(p5,p4,p3,p1,p2, ncol=2)
    glm1$r.squared

    All_lags %>%
      mutate(predicted = predict(glm1, All_lags)$predictions, res = Length - predicted) %>%
      ggplot(aes(predicted, res)) +
      geom_point(colour = "#ff6767", alpha = 0.3) +
      labs(title = "Predicted and observed") +  theme_bw(18) + geom_abline(slope=0)

    All_lags %>%
      mutate(predicted = predict(glm1, All_lags)$predictions, res = Length - predicted) %>%
      ggplot(aes(predicted, Length)) +
      geom_point(colour = "#ff6767", alpha = 0.3) +
      labs(title = "Predicted and observed") +  theme_bw(18) + geom_abline(slope=1)



  	glm2 <- ranger(Latitude ~ ., data=All_lags %>% dplyr::select(! c(Length, julian_release_std)), importance = "impurity")
  	glm2 <- ranger(Latitude ~ ., data=All_lags %>% dplyr::select(! c(Length)), importance = "impurity", mtry=3, num.trees=10000)
  	glm2 <- ranger(Latitude ~ ., data=All_lags %>% dplyr::select(! c(Length)), importance = "impurity", mtry=2, num.trees=10000)
  	glm2 <- ranger(Latitude ~ ., data=All_lags %>% dplyr::select(! c(Length)), importance = "impurity", mtry=1, num.trees=10000)
  	glm2 <- ranger(Latitude ~ ., data=All_lags, importance = "impurity")
  	glm2b <- (lm(Latitude ~ cycle + direction * Release_year, data=All_lags))
  	# glm2b <- glm(Latitude ~ cycle + direction * Release_year, data=All_lags, family = gaussian(link="log"))
  	DHARMa::simulateResiduals(glm2b, plot=TRUE)
  	partial(glm2b, pred.var = "direction", plot=TRUE, rug=TRUE)
  	partial(glm2b, pred.var = "cycle", plot=TRUE, rug=TRUE)
  	partial(glm2b, pred.var = "Release_year", plot=TRUE, rug=TRUE)
  	partial(glm2b, pred.var = c("direction", "Release_year"), plot=TRUE, chull = TRUE)
  	summary(glm2b)
  	anova(glm2b)


  	glm2$prediction.error

  	library(pdp)
  	library(vip)

    p1 <- partial(glm2, pred.var = "direction", plot=TRUE, rug=TRUE)
    p2 <- partial(glm2, pred.var = "cycle", plot=TRUE, rug=TRUE)  # Figure 2
    p3 <- partial(glm2, pred.var = "Release_year", plot=TRUE, rug=TRUE)  # Figure 2
    p4 <- partial(glm2, pred.var = "julian_release_std", plot=TRUE, rug=TRUE, smooth=TRUE)  # Figure 2
    # p4 <- glm2 %>% partial(pred.var = "julian_release_std") %>% autoplot(smooth = TRUE, ylab = "yhat") +
    #   theme_light()
    pp2 <- grid.arrange(p4, p3,p1,p2, ncol=2)
    pp2 <- grid.arrange(p3,p1,p2, ncol=2)
    glm2$r.squared

    p4 <- partial(glm2, pred.var = c("direction", "cycle"), plot=TRUE, chull = TRUE)
    p5 <- partial(glm2, pred.var = c("direction", "Release_year"), plot=TRUE, chull = TRUE)

    All_lags %>%
      mutate(predicted = predict(glm2, All_lags)$predictions) %>%
      ggplot(aes(predicted, Latitude)) +
      geom_point(colour = "#ff6767", alpha = 0.3) +
      labs(title = "Predicted and observed") +  theme_bw(18) + geom_abline(slope=1)
    All_lags %>%
      mutate(predicted = predict(glm2, All_lags)$predictions, res = Latitude - predicted) %>%
      ggplot(aes(predicted, res)) +
      geom_point(colour = "#ff6767", alpha = 0.3) +
      labs(title = "Predicted and observed") +  theme_bw(18) + geom_abline(slope=0)



  # Variance partitioning
  	b <- gam(y ~  1, family=binomial, data = Dec_0lag[[1]]$data)
  	b0 <- Dec_0lag[[1]][[9]]
  	b1 <- gam(y ~  -1 + s(Latitude), sp=b0$sp[1], family=binomial, data = Dec_0lag[[1]]$data)
  	b2 <- gam(y ~  -1 + Catch_year, family=binomial, data = Dec_0lag[[1]]$data)
  	b3 <- gam(y ~  -1 + s(Length), sp=b0$sp[2], family=binomial, data = Dec_0lag[[1]]$data)
  	b4 <- gam(y ~  -1 + s(julian_recapture_std_scaled, k=3), sp=b0$sp[3], family=binomial, data = Dec_0lag[[1]]$data)
    (deviance(b1)-deviance(b0))/deviance(b)
    (deviance(b2)-deviance(b0))/deviance(b)
    (deviance(b3)-deviance(b0))/deviance(b)
    (deviance(b4)-deviance(b0))/deviance(b)

  	with(Feb_0lag[[4]]$data, table(y, Catch_year))
  	b1 <- gam(y ~  -1 + s(Latitude) + Catch_year + s(Length) + factor(month>12), family=binomial, data = Feb_0lag[[4]]$data)
  	b1 <- gam(y ~  -1 + s(Latitude) + Catch_year + s(Length), family=binomial, data = Feb_0lag[[4]]$data)
  	summary(b1)
  	par(mfrow=c(2,2))
  	plot.gam(b1, all.terms=T)





  	# -----------------------------------------------------------------------
  	# ---- Fig 1
  	# -----------------------------------------------------------------------

  	library(ggOceanMaps)
  	library(sf)
  	Norway <- st_read("shapefile/ne_10m_land.shp")
  	Norway <- st_crop(Norway, c(xmin = -34, ymin = 49, xmax = 22, ymax = 72))

  	dt <- data.frame(lon = c(-31, -31, 19, 19), lat = c(45, 70, 70, 45))

  	Catch_data <- taggart::tg_catches(lower=T)
  	Mark_recap_data <- taggart::tg_expeditions(lower=T)
  	d <- Mark_recap_data %>% filter(lubridate::year(releasedate) <2021, lubridate::year(releasedate) >= 2014)
  	d <- filter(d,
  	            ices_rectangle %in% unique(d$ices_rectangle)[unlist(sapply(31:52, function(x) grep(x, unique(d$ices_rectangle))))])

  	d <- filter(d,
  	            !(ices_rectangle %in% unique(d$ices_rectangle)[unlist(sapply("F", function(x) grep(x, unique(d$ices_rectangle))))]))


  	# Releases:
  	df_rel <- d %>% filter(between(lubridate::year(releasedate),2014,2020))
  	df_rel <- mutate(df_rel,
  	                 iceslon =mapplots::ices.rect(ices_rectangle)$lon,
  	                 iceslat = mapplots::ices.rect(ices_rectangle)$lat)


  	# Recaptures:
  	recap_raw <- df_rel %>% filter(!is.na(recapturedate)) %>%
  	  left_join(Catch_data, by = "catchid") %>%
  	  mutate(RecaptureYear = ifelse(lubridate::month(recapturedate)>2,
  	                                lubridate::year(recapturedate),
  	                                lubridate::year(recapturedate)-1),
  	         ReleaseYear = lubridate::year(releasedate),
  	         YearsOut = RecaptureYear - ReleaseYear) %>%
  	  filter(YearsOut<3, between(ReleaseYear,2014,2020),
  	         !is.na(ices_rectangle.y))


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


  	### May current data
  	Current_data <- Dat %>% filter(time %in% c(564276, 573036, 581820, 590580, 599340, 608100, 616884))
  	Current_data$Release_year <- (2014:2020)[match(Current_data$time, c(564276, 573036, 581820, 590580, 599340, 608100, 616884))]
  	Current_data <- Current_data %>% filter(lat>51, lon<15, lat<68.5)

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

  	new_sf <- ICESecoregion %>% filter(Ecoregion %in% c("Norwegian Sea", "Icelandic Waters", "Greater North Sea", "Celtic Seas")) %>%
  	  group_by(Ecoregion) %>% sfheaders::sf_remove_holes()



  	g1 <- ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
  	  geom_sf(data = Bathy, aes(col=depth, fill=depth)) +
  	  scale_colour_manual(values=cc, name = "Depth") + scale_fill_manual(values=cc, name = "Depth") +
  	  geom_sf(data=Norway, col=grey(0.3)) +
  	  # new_scale_color() +
  	  geom_sf(data=new_sf, aes(geometry = geometry), col = "black", fill=NA, size=0.9) +
  	  geom_text(data = data.frame(x=c(-19,0,2.3,-12), y=c(62,68,56.5,55), label=c("Icelandic \nWaters", "Norwegian \nSea", "Greater \nNorth \nSea", "Celtic\nSeas")),
  	            aes(x=x, y=y, label=label), col="white", size=5, fontface=2) +
  	  geom_point(data = df_rel %>% distinct(longitude,latitude),    aes(x = longitude, y = latitude), col = "red", alpha = 0.5, size = .5)+
  	  geom_point(data = recap_raw %>% distinct(clon,clat), aes(x = clon, y = clat), col = "yellow", size = .6)+
  	  coord_sf(xlim=c(-31, 19), ylim=c(50, 70)) +
  	  labs(x="Longitude", y="Latitude") +
  	  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
  	                     legend.title = element_text(face = "bold", size = 12),
  	                     strip.background = element_rect(fill = "white"),
  	                     strip.text = element_text(face ="bold", size = 14),
  	                     axis.text = element_text(size =13),
  	                     axis.title = element_text(size =14))
  	#panel.spacing.x = unit(6.5, "mm"))+
  	# ggsave(g1, filename = "MS/figs/Fig2_May_Current_bathym.pdf",
  	#        width=28, height=32, units="cm", dpi = 450)
  	#


  	ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
    	  geom_sf(data=bla, aes(geometry = geometry, col=Ecoregion), fill=NA, size=2)


  	g2 <- ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
  	  geom_sf(data = Bathy, aes(col=depth, fill=depth), show.legend = FALSE) +
  	  scale_colour_manual(values=cc, name = "Depth") + scale_fill_manual(values=cc, name = "Depth") + geom_sf(data=Norway, col=grey(0.3)) +
  	  coord_sf(xlim=c(-31, 19), ylim=c(50, 70)) +
  	  new_scale_color() +
  	  new_scale_fill() +
  	  geom_segment(data = Current_data %>% filter(speed > 0.15, Release_year %in% c(2014, 2018)), aes(x=X, y=Y, xend = X + EW_comp*5, yend = Y + NS_comp*5, col=speed),
  	               arrow = arrow(angle = 15, length = unit(0.03, "inches"), type = "closed"), alpha = 0.8, size=0.3) +
  	  scale_color_gradient(low="lightpink1", high="darkred", name = "Current\nspeed (m/s)") +
  	  facet_wrap(~ Release_year, ncol=2) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
  	                     legend.title = element_text(face = "bold", size = 12),
  	                     strip.background = element_rect(fill = "white"),
  	                     strip.text = element_text(face ="bold", size = 14),
  	                     axis.text = element_text(size =13),
  	                     axis.title = element_text(size =14)) +
  	  labs(x="Longitude", y="Latitude")


  pp <- gridExtra::grid.arrange(g1, g2, nrow=2)
  ggsave(pp, file="MS/Figs/Fig1.png", width=10, height=10, dpi=400)



