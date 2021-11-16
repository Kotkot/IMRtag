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

Norway <- st_read("C:/Users/a23092/Dropbox/IMR_projects/Shapefiles/ne_10m_land.shp")
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

	model_type <- c("Iceland", "Norway", "Northsea")[2]
	run_directionality <- function(model_type, month, best="AIC", data_origin="0year"){
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
  			         direction = ifelse(cLon < -10 & cLat>=62, 1L, ifelse((cLon >= -10 & cLat>=62), 2L, ifelse((cLon >= -10 & cLat<62), 0L, 3L)))
  			         )

  			Data_mackerel_use_Ireland_select %>% group_by(toiceland, Catch_year) %>% summarize(n=n())

  			Data_mackerel_use_Ireland_select %>% filter(toiceland == 1) %>% ggplot(aes(x=Latitude)) + geom_histogram()


  			hist(Data_mackerel_use_Ireland_select$dist_toiceland, breaks=50, xlab="Distance from Norway coast (m)", main="")

  	## Analysis of the EW movement
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
  			# m0 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled),
  			#                ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  # 			#           family=multinom(K=2), data = www)  # failed
  # 			m1 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			m2 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			# m3 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled),
  # 			#                ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  # 			#           family=multinom(K=2), data = www) # failed
  # 			m4 <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			m5 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			m6 <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			m7 <- gam(list(direction ~  -1 + s(Latitude) + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			m8 <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year, bs="ts") + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  # 			mall <- gam(list(direction ~  -1 + s(Latitude, by=Catch_year, bs="ts") + Catch_year + s(Length) + s(julian_recapture_std_scaled),
  # 			               ~  -1 + s(Latitude, by=Catch_year, bs="ts") + Catch_year + s(Length, by=Catch_year, bs="ts") + s(julian_recapture_std_scaled)),
  # 			          family=multinom(K=2), data = www)
  #       AICs <- AIC(m1,m2,m4,m5,m6,m7,m8,mall)
  #       BICs <- BIC(m1,m2,m4,m5,m6,m7,m8)
  #
  #       gam.selection()
  #       if (best == "AIC") best_model <- list(m1,m2,m4,m5,m6,m7,m8)[[which.min(AICs$AIC)]]
  #       if (best == "Manual") best_model <- m6
  #       plot.gam(mall, residuals=TRUE, pch=19, cex=0.5, scale=0,pages = 2)
  #       plot.gam(best_model, pch=19, cex=0.5, scale=0,pages = 1)
  #       # qq.gam(best_model, rep = 100, level = 0.9, type = "deviance", rl.col = 2,  rep.col = "gray80", cex=3)
  #       pred <- predict(best_model, type="response")
  #       # p1 <- qq.gam(best_model, rep = 5000, level = 0.9, type = "deviance", rl.col = 2,  rep.col = "gray80", cex=3, s.rep=5000)
  #       # points(qnorm(as.numeric(str_remove(names(p1), pattern= "%"))/100), as.numeric(p1), col="red")
  #       # DHARMa::simulateResiduals(best_model, plot=TRUE)
  #       u1 <- sapply(which(www$toiceland==1), function(x) pbinom(q=www$toiceland[x], size=1, prob=pred[x,2]))
  #       u2 <- pbinom(q=www$to_norway[which(www$to_norway==1)], size=1, prob=pred[which(www$to_norway==1),3])
  #       u3 <- pbinom(q=www$to_northsea[which(www$to_northsea==1)], size=1, prob=pred[which(www$to_northsea==1),1])
  #       u <- c(u1,u2,u3)
  #
  #       u <- pbinom(q=www$toiceland, size=1, prob=pred[,2])
  #       u <- pbinom(q=www$to_norway, size=1, prob=pred[,3])
  #       u <- pbinom(q=www$to_northsea, size=1, prob=pred[,1])
  #       hist(u)
  #       qqnorm(qq)
  #       qq <- qnorm(u)
  #       range(qq)
  #
  #       table(apply(pred, 1, which.max), www$direction+1)
  #       aa <- table(apply(predict(m1, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #       aa <- table(apply(predict(m2, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #       aa <- table(apply(predict(m4, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #       aa <- table(apply(predict(m5, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #       aa <- table(apply(predict(m6, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #       aa <- table(apply(predict(m7, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #       aa <- table(apply(predict(m8, type="response"), 1, which.max), www$direction+1); diag(aa) <- 0; sum(aa)
  #

  			### 10-fold cross-validation

  			library(groupdata2)
  			library(furrr)
  			library(progressr)

  			kfolding <- function(it, kfolds = 10, model_type="Iceland"){

  			  set.seed(it)
  			  ERRORS <- c()
  			  AICs <- c()
  			  BICs <- c()
  			  if(model_type == "Iceland") www$y <- www$toiceland
  			  if(model_type == "Norway") www$y <- www$to_norway
  			  if(model_type == "Northsea") www$y <- www$to_northsea
  			  www1 <- www %>% ungroup() %>% fold(k = kfolds, method="n_rand")

  			  for (k in 1:kfolds){

  			    testing <- www1[www1$.folds == k,]
  			    training <- www1[-testing$ID,]

  			    if (kfolds == 1) dat_use <- www
  			    if (kfolds > 1) dat_use <- training

  			    m0=m1=m2=m3=m4=NULL
  			    m0 <- (gam(y ~  1, family=binomial,
  			              data = dat_use))
  			    m1 <- (gam(y ~  -1 + s(Latitude) + Catch_year  + s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
  			              data = dat_use))
  			    m2 <- (gam(y ~  -1 + s(Latitude) + Catch_year + s(Length, by=Catch_year) + s(julian_recapture_std_scaled, k=3), family=binomial,
  			              data = dat_use))
  			    m3 <- (gam(y ~  -1 + s(Latitude, by=Catch_year) + Catch_year + s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
  			              data = dat_use))
  			    m4 <- (gam(y ~  -1 + s(Latitude) + Catch_year  + s(Length) + s(julian_recapture_std_scaled, by=Catch_year), family=binomial,
  			              data = dat_use))
  			    mse0 = mse1 = mse2 = mse3 = mse4 = NA

  			    if (is(m0) != "try-error")   pred0 <- predict(m0, type="response", newdata = testing)
  			    if (is(m1) != "try-error")   pred1 <- predict(m1, type="response", newdata = testing)
  			    if (is(m2) != "try-error")   pred2 <- predict(m2, type="response", newdata = testing)
  			    if (is(m3) != "try-error")   pred3 <- predict(m3, type="response", newdata = testing)
  			    if (is(m4) != "try-error")   pred4 <- predict(m4, type="response", newdata = testing)

  			    tonum <- function(x) as.numeric(as.character(x))
  			    obs <- tonum(testing$y)
  			    if (is(m0) != "try-error")   mse0 <- mean(apply(pred0-obs, 1, function(x) sum(x^2)))
  			    if (is(m1) != "try-error")   mse1 <- mean(apply(pred1-obs, 1, function(x) sum(x^2)))
  			    if (is(m2) != "try-error")   mse2 <- mean(apply(pred2-obs, 1, function(x) sum(x^2)))
  			    if (is(m3) != "try-error")   mse3 <- mean(apply(pred3-obs, 1, function(x) sum(x^2)))
  			    if (is(m4) != "try-error")   mse4 <- mean(apply(pred4-obs, 1, function(x) sum(x^2)))

  			    nAICs <- AIC(m0,m1,m2,m3,m4)
  			    nBICs <- BIC(m0,m1,m2,m3,m4)
  			    error <- c(mse0, mse1, mse2, mse3, mse4)
  			    ERRORS <- rbind(ERRORS, error)
  			    AICs <- rbind(AICs, nAICs)
  			    BICs <- rbind(BICs, nBICs)
  			  }

  			  RESULT <- list(m0=m0,m1=m1,m2=m2,m3=m3,m4=m4,ERRORS=ERRORS, AICs=AICs, BICs=BICs)

  			  if(kfolds == 1) return(RESULT)
  			  if(kfolds > 1) return(ERRORS)
  			}

  			### Go to Iceland

        if (model_type == "Iceland"){

          main_iceland = kfolding(it=1, kfolds=1, model_type="Iceland")

          CVs <- kfolding(it=1, kfolds = 10, model_type="Iceland")
          apply(CVs, 2, mean, na.rm=T)
          CVs <- kfolding(it=20, kfolds = 10, model_type="Iceland")
          apply(CVs, 2, mean, na.rm=T)

          main_iceland$AICs
          concurvity(main_iceland$m1)
          concurvity(main_iceland$m2)
          concurvity(main_iceland$m3)
          concurvity(main_iceland$m4)

          library(DHARMa)
          library(ggeffects)
          par(mfrow=c(2,2))
          plot.gam(main_iceland$m1, all.terms=T)
          simul <- simulateResiduals(main_iceland$m1, plot=TRUE)
          par(mfrow=c(2,2))
          plot(www$Latitude, simul$scaledResiduals)
          plot(www$Length, simul$scaledResiduals)
          plot(www$julian_recapture_std, simul$scaledResiduals)
          plot(www$Catch_year, simul$scaledResiduals)

          mydf <- ggpredict(main_iceland$m1, terms = "Latitude", back.transform = FALSE)
          ggplot(mydf, aes(x, predicted)) +
            geom_line() +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
  			}

  			if (model_type == "Norway"){
  			  main_norway = kfolding(it=1, kfolds=1, model_type="Norway")
  			  main_norway$AICs
  			  main_norway$BICs
  			  concurvity(main_norway$m1)
  			  concurvity(main_norway$m2)
  			  concurvity(main_norway$m3)
  			  concurvity(main_norway$m4)
  			  library(DHARMa)
  			  simul <- simulateResiduals(main_norway$m1, plot=TRUE)
  			  par(mfrow=c(2,2))
  			  plot.gam(main_norway$m1, all.terms=TRUE)

  			  CVs <- kfolding(it=1, kfolds = 10, model_type="Norway")
  			  apply(CVs, 2, mean, na.rm=T)
  			  CVs <- kfolding(it=20, kfolds = 10, model_type="Norway")
  			  apply(CVs, 2, mean, na.rm=T)

          par(mfrow=c(2,2))
  			  plot(www$Latitude, simul$scaledResiduals)
  			  plot(www$Length, simul$scaledResiduals)
  			  plot(www$julian_recapture_std, simul$scaledResiduals)
  			  plot(www$Catch_year, simul$scaledResiduals)

  			}

  			if (model_type == "Northsea"){
  			  main_northsea = kfolding(it=1, kfolds=1, model_type="Northsea")
  			  main_northsea$AICs
  			  main_northsea$BICs
  			  concurvity(main_northsea$m1)
  			  concurvity(main_northsea$m2)
  			  concurvity(main_northsea$m3)
  			  concurvity(main_northsea$m4)
  			  library(DHARMa)
  			  simul <- simulateResiduals(main_northsea$m1, plot=TRUE)
  			  par(mfrow=c(2,2))
  			  plot.gam(main_northsea$m1, all.terms=TRUE)

  			  CVs <- kfolding(it=1, kfolds = 10, model_type="Northsea")
  			  apply(CVs, 2, mean, na.rm=T)
  			  CVs <- kfolding(it=20, kfolds = 10, model_type="Northsea")
  			  apply(CVs, 2, mean, na.rm=T)

  			  library(DHARMa)
  			  par(mfrow=c(2,2))
  			  simul <- simulateResiduals(main_northsea$m1, plot=TRUE)
  			  par(mfrow=c(2,2))
  			  plot(www$Latitude, simul$scaledResiduals)
  			  plot(www$Length, simul$scaledResiduals)
  			  plot(www$julian_recapture_std, simul$scaledResiduals)
  			  plot(www$Catch_year, simul$scaledResiduals)
  			}
  			AICs <- AIC(m1,m2,m4,m5,m6)

  			if (best == "AIC") best_model <- list(m1,m2,m4,m5,m6)[[which(AICs == min(AICs))]]
  			if (best == "Manual") best_model <- m1
  			summary(best_model)
  			plot.gam(best_model, residuals=TRUE, pch=19, cex=0.5, scale=0,pages = 1)
  			plot.gam(best_model, pch=19, cex=0.5, scale=0,pages = 1)
  			# qq.gam(m1, rep = 0, level = 0.9, type = "deviance", rl.col = 2,  rep.col = "gray80", cex=3, s.rep=5000)
  			# p1 <- qq.gam(best_model, rep = 5000, level = 0.9, type = "deviance", rl.col = 2,  rep.col = "gray80", cex=3, s.rep=5000)
  			# points(qnorm(as.numeric(str_remove(names(p1), pattern= "%"))/100), as.numeric(p1), col="red")
  			plot(DHARMa::simulateResiduals(best_model),pch=".")
  			dev.copy2pdf(file=paste0(getwd(), "/MS/figs/EW_mvt/resid_", label, model_type, data_origin, ".pdf"), out.type="pdf")
  			dev.off()


  			www$res <- residuals(best_model, type="deviance")
  			p1 <- ggplot(www, aes(x=Length, y=res)) + geom_point() + geom_smooth() + theme_bw()
  			p2 <- ggplot(www, aes(x=Latitude, y=res)) + geom_point() + geom_smooth() + theme_bw()
  			p3 <- ggplot(www, aes(x=julian_release_std, y=res)) + geom_point() + geom_smooth() + theme_bw()
  			p4 <- ggplot(www, aes(x=julian_recapture_std_scaled, y=res)) + geom_point() + geom_smooth() + theme_bw()
  			p5 <- ggplot(www, aes(x=factor(Catch_year), y=res)) + geom_boxplot() + theme_bw()
  			p6 <- ggplot(www, aes(x=factor(Release_timing_fact), y=res)) + geom_boxplot() + theme_bw()

  			p0 <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow=3)
  			ggsave(p0, file=paste0(getwd(), "/MS/figs/EW_mvt/residuals",label, model_type, data_origin, ".pdf"),width=20, height=14, units="cm", dpi = 450)



  			# nsamp <- www %>% group_by(Release_timing_fact, toiceland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			#   ungroup() %>% group_by(Release_timing_fact, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			#                                            ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			#   filter(toiceland_bin == 1)
  			if (model_type == "Iceland"){
  			  nsamp <- www %>% group_by(toiceland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toiceland_bin == 1)
  			}
  			if (model_type == "Norway"){
  			    nsamp <- www %>% group_by(tonorway_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorway_bin == 1)
  			 }
  			if (model_type == "Northsea"){
  			    nsamp <- www %>% group_by(tonorthsea_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  mutate(n_tot=sum(n),
  			         ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorthsea_bin == 1)
  			 }

  			pp <- visreg::visreg(best_model, "Latitude", plot=FALSE)
  			p0 <- ggplot(pp$fit, aes(x=Latitude, y=visregFit)) + geom_line() +
  			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
  			  # facet_wrap(~Release_timing_fact, scales="free_y") + theme_bw() +
  			   theme_bw() +
  			  ylab("Marginal effect") +
  			 xlab("Latitude") +
  			  theme(axis.title = element_text(size=13),
  			        axis.text = element_text(size=10),
  			        strip.text.x = element_text(size = 10))
  			if (model_type == "Iceland") pp0 <- p0 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==0), aes(x=Latitude)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==1), aes(x=Latitude), col="red") +
  			  geom_text(data=nsamp, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			if (model_type == "Norway") pp0 <- p0 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==0), aes(x=Latitude)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==1), aes(x=Latitude), col="red") +
  			  geom_text(data=nsamp, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			if (model_type == "Northsea") pp0 <- p0 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=Latitude)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=Latitude), col="red") +
  			  geom_text(data=nsamp, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)

  			# ggsave(pp0, file=paste0(getwd(), "/MS/figs/EW_mvt/Marginal_latitude_limit",label, ".pdf"),width=20, height=14, units="cm", dpi = 450)

  			pp <- visreg::visreg(best_model, "julian_recapture_std_scaled", plot=FALSE)
  			p1 <- ggplot(pp$fit, aes(x=julian_recapture_std_scaled, y=visregFit)) + geom_line() +
  			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
  			  # facet_wrap(~Release_timing_fact, scales="free_y") + theme_bw() +
  			   theme_bw() +
  			  ylab("") +
  			  xlab("Recapture date") +
  			  theme(axis.title = element_text(size=13),
  			        axis.text = element_text(size=10),
  			        strip.text.x = element_text(size = 10))
  			if (model_type == "Iceland") pp1 <- p1 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==0), aes(x=julian_recapture_std_scaled)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==1), aes(x=julian_recapture_std_scaled), col="red") +
  			  geom_text(data=nsamp, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			if (model_type == "Norway") pp1 <- p1 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==0), aes(x=julian_recapture_std_scaled)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==1), aes(x=julian_recapture_std_scaled), col="red") +
  			  geom_text(data=nsamp, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			if (model_type == "Northsea") pp1 <- p1 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=julian_recapture_std_scaled)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=julian_recapture_std_scaled), col="red") +
  			  geom_text(data=nsamp, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)

  			# ggsave(pp1, file=paste0(getwd(), "/MS/figs/EW_mvt/Marginal_latitude_limit",label, ".pdf"),width=20, height=14, units="cm", dpi = 450)

  			if (model_type == "Iceland"){
  			  nsamp_length <- www %>% group_by(Catch_year, toiceland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Catch_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                        ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toiceland_bin == 1)
  			}
  			if (model_type == "Norway"){
  			  nsamp_length <- www %>% group_by(Catch_year, tonorway_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Catch_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                        ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorway_bin == 1)
  			}
  			if (model_type == "Northsea"){
  			  nsamp_length <- www %>% group_by(Catch_year, tonorthsea_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Catch_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                        ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorthsea_bin == 1)
  			}

  			pp <- visreg::visreg(best_model, "Length", by ="Catch_year", plot=FALSE)
  			p2 <- ggplot(pp$fit, aes(x=Length, y=visregFit)) + geom_line() +
  			  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
  			  facet_wrap(~Catch_year, scales="free_y") + theme_bw() +
  			  ylab("") +
  			  xlab("Length (cm)") +
  			  theme(axis.title = element_text(size=13),
  			        axis.text = element_text(size=10),
  			        strip.text.x = element_text(size = 10))
  			if (model_type == "Iceland") pp2 <- p2 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==0), aes(x=Length)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(toiceland_bin==1), aes(x=Length), col="red") +
  			  geom_text(data=nsamp_length, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			if (model_type == "Norway") pp2 <- p2 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==0), aes(x=Length)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorway_bin==1), aes(x=Length), col="red") +
  			  geom_text(data=nsamp_length, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  			if (model_type == "Northsea") pp2 <- p2 + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=Length)) +
  			  geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=Length), col="red") +
  			  geom_text(data=nsamp_length, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)

  			if (model_type == "Iceland"){
  			  nsamp_year <- www %>% group_by(Catch_year, toiceland_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Catch_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(toiceland_bin == 1)
  			}
  			if (model_type == "Norway"){
  			  nsamp_year <- www %>% group_by(Catch_year, tonorway_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Catch_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorway_bin == 1)
  			}
  			if (model_type == "Northsea"){
  			  nsamp_year <- www %>% group_by(Catch_year, tonorthsea_bin, .drop=FALSE) %>% summarize(n=n()) %>%
  			  ungroup() %>% group_by(Catch_year, .drop=FALSE) %>% mutate(n_tot=sum(n),
  			                                                             ntot = paste0("n=", n_tot, "(", n, ")")) %>%
  			  filter(tonorthsea_bin == 1)
  			}

  	    pp <- visreg::visreg(best_model, "Catch_year", plot=FALSE)
  			p3 <- ggplot(pp$fit, aes(x=Catch_year, y=visregFit)) + geom_point() +
  			  geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col=grey(0.5), width=0.4) +
  			  theme_bw() +
  			  ylab("") +
  			  xlab("Year") +
  			  theme(axis.title = element_text(size=13),
  			        axis.text = element_text(size=10),
  			        strip.text.x = element_text(size = 10))

        ppp <- grid.arrange(pp0,pp1,pp2,p3,layout_matrix=matrix(c(1,2,3,4,4,4), byrow=T, ncol=3, nrow=2), widths=c(1,1,2), heights=c(1,0.5))
        ggsave(ppp, file=paste0(getwd(), "/MS/figs/EW_mvt/Marginal_effects",label, model_type, data_origin, ".pdf"),width=20, height=20, units="cm", dpi = 450)

        return(ppp)
  	  }

  # same year recapture
  	test1 <- run_directionality("Iceland", month=12, best="AIC", data_origin="0year")
  	test2 <- run_directionality("Norway", month=12, best="AIC", data_origin="0year")
  	test3 <- run_directionality("Northsea", month=12, best="AIC", data_origin="0year")
    pall <- grid.arrange(test1,test2,test3, nrow=2)
    ggsave(pall, file=paste0(getwd(), "/MS/figs/EW_mvt/all_year0_month12.pdf"),width=35, height=29, units="cm", dpi = 450)
  # 1 year after recapture
  	test1 <- run_directionality("Iceland", month=12, best="AIC", data_origin="1year")
  	test2 <- run_directionality("Norway", month=12, best="AIC", data_origin="1year")
  	test3 <- run_directionality("Northsea", month=12, best="AIC", data_origin="1year")
    pall <- grid.arrange(test1,test2,test3, nrow=2)
    ggsave(pall, file=paste0(getwd(), "/MS/figs/EW_mvt/all_year1_month12.pdf"),width=35, height=29, units="cm", dpi = 450)
  # 2 year after recapture
  	test1 <- run_directionality("Iceland", month=12, best="AIC", data_origin="2year")
  	test2 <- run_directionality("Norway", month=12, best="AIC", data_origin="2year")
  	test3 <- run_directionality("Northsea", month=12, best="AIC", data_origin="2year")
    pall <- grid.arrange(test1,test2,test3, nrow=2)
    ggsave(pall, file=paste0(getwd(), "/MS/figs/EW_mvt/all_year2_month12.pdf"),width=35, height=29, units="cm", dpi = 450)

	## This makes me think that I might need to develop a changepoint model (with K components)
			## I sometimes call it "mixture" model below but it is not a mixture model but a changepoint model
			## A Bayesian change point model has therefore been developped below

			# Choice of the design matrix
			m1 <- lm(log_rate ~ factor(Tag_area)*factor(Release_timing) + factor(Release_year) + length_scaled + julian_recapture_scaled, data=Data_mackerel_use_Ireland_select)
			m_frame <- model.frame(m1)
			XX <- model.matrix(m1, m_frame)

			# Choice of threshold values: Either the travel distance or the time of the year
			yyy <- Data_mackerel_use_Ireland_select$log_rate
			Data_mackerel_use_Ireland_select$julian_release <-  as.numeric(julian(Data_mackerel_use_Ireland_select$ReleaseDate, as.POSIXct(paste0(2014, "-01-01"), tz = "GMT")))
			Data_mackerel_use_Ireland_select$julian <-  as.numeric(julian(Data_mackerel_use_Ireland_select$RecaptureDate, as.POSIXct(paste0(2014, "-01-01"), tz = "GMT")))
			Data_mackerel_use_Ireland_select$julian_std <-  Data_mackerel_use_Ireland_select$julian %% 365
			yyy1 <- Data_mackerel_use_Ireland_select$julian_std
			threshold_vals <- as.numeric(quantile(yyy1, seq(0.1, 0.9, by=0.05)))
			threshold_vals <- sort(yyy1[which(yyy1>quantile(yyy1,0.1) & yyy1<quantile(yyy1,0.9))])
			threshold_vals_group <- cut(Data_mackerel_use_Ireland_select$julian_std, c(0, threshold_vals, 365))
			threshold_vals_group <- as.numeric(as.character(factor(threshold_vals_group, labels=1:(length(threshold_vals)+1))))
			threshold_vals_group_start <- c(1, which(diff(threshold_vals_group, lag=1) == 1)+1)
			threshold_vals_group_end <- c(which(diff(threshold_vals_group, lag=1) == 1), length(threshold_vals_group))


      if(do_bayesian == TRUE){

        # With 2 groups
        Data <- list(
          K=2,  # number of mixture components
          N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
          Nx=ncol(XX),   # the fixed effect part
          X=XX,          # the design matrix for the fixed effect
          Nthres=length(threshold_vals),
          thresh=threshold_vals,
          thresh_start=threshold_vals_group_start,
          thresh_end=threshold_vals_group_end,
          mean_diff_tag_area= mean_diff_tag_area,
          is_from_South=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "South_Ireland",1,0),
          y = Data_mackerel_use_Ireland_select$log_rate
        )

        library(rstan)
        options(mc.cores = 3)
        rstan_options(auto_write = TRUE)
        mixture_mod <- stan( file = paste0(getwd(), "/src/mackerel_mvt_model_threoshold.stan"), data = Data,
                             iter = 20000 , warmup=15000 , chains = 3, thin=10,
                             control = list(adapt_delta = 0.99, max_treedepth = 15), seed=123)
        # the use of dynamic programming to linearise the problem... does not work too well...
        mixture_moda <- stan( file = paste0(getwd(), "/src/mackerel_mvt_model_threoshold_dp.stan"), data = Data,
                              iter = 5000 , warmup=3500 , chains = 3, thin=10,
                              control = list(adapt_delta = 0.99, max_treedepth = 20), seed=123)

        # With 3 groups (computer is struggling)
        Data3 <- list(
          K=3,  # number of mixture components
          N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
          Nx=ncol(XX),   # the fixed effect part
          X=XX,          # the design matrix for the fixed effect
          Nthres=length(threshold_vals),
          thresh=threshold_vals,
          thresh_start=threshold_vals_group_start,
          thresh_end=threshold_vals_group_end,
          mean_diff_tag_area= mean_diff_tag_area,
          is_from_South=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "South_Ireland",1,0),
          y = Data_mackerel_use_Ireland_select$log_rate
        )

        set.seed(1)
        inits1 <- list(beta=matrix(c(rep(10,3),runif(9,-2,2),runif(6*3,-0.05,0.05)),byrow=T, ncol=3), sigma=runif(1,0.1,0.5))
        inits2 <- list(beta=matrix(c(rep(10,3),runif(9,-2,2),runif(6*3,-0.05,0.05)),byrow=T, ncol=3), sigma=runif(1,0.1,0.5))
        inits3 <- list(beta=matrix(c(rep(10,3),runif(9,-2,2),runif(6*3,-0.05,0.05)),byrow=T, ncol=3), sigma=runif(1,0.1,0.5))

        inits_gen <- function(chain_id = 1) {
          list(beta=matrix(c(rep(10,3),runif(9,-2,2),runif(6*3,-0.05,0.05)),byrow=T, ncol=3),
               sigma=runif(3,0.1,0.5)) #,
          # mu=XX %*% inits1$beta,
          # lp=rep(0, length(threshold_vals)))
        }
        init_ll <- lapply(1, function(id) inits_gen(chain_id = id))


        library(rstan)
        options(mc.cores = 3)
        rstan_options(auto_write = TRUE)
        mixture_mod3 <- stan( file = paste0(getwd(), "/src/mackerel_mvt_model_threshold_2breaks.stan"), data = Data3,
                              iter = 5000 , warmup=3500 , chains = 3, thin=10,
                              control = list(adapt_delta = 0.99, max_treedepth = 15), seed=123)
        # the use of dynamic programming to linearise the problem... does not work too well...
        mixture_mod3b <- stan( file = paste0(getwd(), "/src/mackerel_mvt_model_threshold_2breaks_dp.stan"), data = Data3,
                               iter = 500, warmup=350, chains = 1, thin=10,
                               #init=init_ll,
                               control = list(adapt_delta = 0.99, max_treedepth = 20), seed=123)

        # Some model checking + residual analysis of the Bayesian model
        library(shinystan)
        launch_shinystan(mixture_mod)
        check_hmc_diagnostics(mixture_mod)     # need no warning
        get_low_bfmi_chains(mixture_mod)
        stan_diag(mixture_mod)
        stan_rhat(mixture_mod)

        sampler_params <- get_sampler_params(mixture_mod, inc_warmup = FALSE)
        sampler_params_chain1 <- sampler_params[[1]]
        colnames(sampler_params_chain1)
        mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
        print(mean_accept_stat_by_chain)
        max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
        print(max_treedepth_by_chain)

        summary(mixture_mod)$summary

        LP <- extract(mixture_mod, pars="lp")
        LP_mean <- apply(exp(LP$lp), 2, mean)

        windows()
        plot(Data$thresh, LP_mean)
        best_LP <- which(LP_mean == max(LP_mean))
        (threshold <- Data$thresh[best_LP])

        par_list <- c("y_gen")
        y_pred <- rstan::extract(mixture_mod, pars=par_list,inc_warmup=FALSE, permuted=FALSE)
        #essai <- array(y_pred[,1,], dim=c(200,17,903))[,best_LP,]
        y_pred_new <- rbind(apply(array(y_pred[,1,], dim=c(200,17,903)), c(1,3), function(x) sum(x*LP_mean/sum(LP_mean))),
                            apply(array(y_pred[,2,], dim=c(200,17,903)), c(1,3), function(x) sum(x*LP_mean/sum(LP_mean))),
                            apply(array(y_pred[,3,], dim=c(200,17,903)), c(1,3), function(x) sum(x*LP_mean/sum(LP_mean))))
        y_pred_summary <- data.frame(n=rep(seq(1,250),4)[1:length(Data$y)], ID = seq_along(Data$y), obs = Data$y,mean=apply(y_pred_new,2,mean), l95=apply(y_pred_new, 2, function(x) quantile(x,0.025)), u95=apply(y_pred_new, 2, function(x) quantile(x,0.975)))
        y_pred_summary$Cat <- cut(y_pred_summary$ID, seq(0,1000,by=250))
        head(y_pred_summary)
        ggplot(y_pred_summary, aes(x=n, y=obs)) + geom_point(col="red") +
          geom_errorbar(aes(x=n, ymin=l95, ymax=u95)) + theme_bw() + facet_grid(Cat~., scales="free")

        # residual analysis
        plot(y_pred_summary$ID, (y_pred_summary$mean-y_pred_summary$obs)/y_pred_summary$mean)
        abline(h=0, lty=2)
        qqnorm(y=(y_pred_summary$mean-y_pred_summary$obs)/sd((y_pred_summary$mean-y_pred_summary$obs)))
        abline(0,1, lty=2)

        # Plot marginal effects
        Beta_estimates <- summary(mixture_mod)$summary[grep("beta", rownames(summary(mixture_mod)$summary)),]
        group1 <- Beta_estimates[seq(1,20,by=2),]
        group2 <- Beta_estimates[seq(2,20,by=2),]
        aaa <- cbind(group1[,1], group2[,1])
        rownames(aaa) <- colnames(XX)
        print(aaa)



      }


		## Now doing the same model but using TMB

  			# Some plotting to justify covariate parametrization
  			par(mfrow=c(3,2), oma=c(1,3,1,1), mar=c(4,2,1,1), cex.lab=1.2, cex.axis=1.2, cex.main=1.3)
  			with(Data_mackerel_use_Ireland_select, boxplot(log_rate ~ Release_year, xlab="Release year", ylab="log(movement rate)"))
  			with(Data_mackerel_use_Ireland_select, boxplot(log_rate ~ Length, xlab="Fish body length (cm)", ylab="log(movement rate)"))
  			with(Data_mackerel_use_Ireland_select, plot(julian_recapture_std, log_rate, xlab="Recapture date (julian days)", ylab="log(movement rate)", ylim=c(8, 11.6)))
  			abline(v=212, lty=2)
  			abline(v=243, lty=2)
  			abline(v=273, lty=2)
  			abline(v=304, lty=2)
  			text(x=198, y=11, "July", srt=45)
  			text(x=228, y=11, "August", srt=45)
  			text(x=258, y=11, "September", srt=45)
  			text(x=288, y=11, "October", srt=45)
  			text(x=318, y=11, "November", srt=45)
  			with(Data_mackerel_use_Ireland_select, boxplot(log_rate ~ Release_timing, xlab="Release timing", ylab="log(movement rate)", names=c("<May22nd", "≥May22nd")))
  			#with(Data_mackerel_use_Ireland_select, boxplot(log_rate ~ as.character(Tag_area), xlab="Latitude", ylab="log(movement rate)"))
  			with(Data_mackerel_use_Ireland_select, plot(log_rate ~ Latitude, xlab="Latitude (release)", ylab="log(movement rate)"))
  			mtext(side=2, line=1.5, "log(linear movement rate (m/day))", outer=T)
  			grid::grid.text("a)", x=unit(0.05, "npc"), y=unit(0.97, "npc"))# with(Data_mackerel_use_Ireland_select, plot(length, log_rate))
  			grid::grid.text("b)", x=unit(0.52, "npc"), y=unit(0.97, "npc"))# with(Data_mackerel_use_Ireland_select, plot(length, log_rate))
  			grid::grid.text("c)", x=unit(0.05, "npc"), y=unit(0.66, "npc"))# with(Data_mackerel_use_Ireland_select, plot(length, log_rate))
  			grid::grid.text("d)", x=unit(0.52, "npc"), y=unit(0.66, "npc"))# with(Data_mackerel_use_Ireland_select, plot(length, log_rate))
  			grid::grid.text("e)", x=unit(0.05, "npc"), y=unit(0.35, "npc"))# with(Data_mackerel_use_Ireland_select, plot(length, log_rate))

  			with(Data_mackerel_use_Ireland_select, plot(log_rate ~ Latitude, xlab="Release point latitude", ylab="Log movement rate"))
  			with(subset(Data_mackerel_use_Ireland_select, Release_timing=="First_half"), plot(log_rate ~ Latitude))
  			with(subset(Data_mackerel_use_Ireland_select, Release_timing=="Second_half"), points(log_rate ~ Latitude, col="red"))

  			ggplot(subset(Data_mackerel_use_Ireland_select, Release_year=2014), aes(x=Latitude, y=log_rate, col=julian_std)) + geom_point()

  			# map_area <- get_map(location= c(left = -35, bottom = 51, right = 21, top = 72))
  			# map_area <- get_map(location= c(left = -15, bottom = 51, right = -8, top = 55))
  			ggmap(map_area) + geom_jitter(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat, col=julian_recapture_std, pch=Tag_area)) +
  			 scale_color_viridis_c()

  			ggmap(map_area) + geom_jitter(data=Data_mackerel_use_Ireland_select, aes(x=Longitude, y=Latitude), cex=0.1)

  			ggplot(data=Data_mackerel_use_Ireland_select %>% filter(julian_recapture_std < 270)) +
  			         geom_sf(data = Norway) + theme_bw() +
  			  geom_jitter(aes(x=cLon, y=cLat, col=Length),size=1) +
  			  scale_color_viridis_c(option="C")

  			ggplot(data=Data_mackerel_use_Ireland_select %>%
  			         filter(julian_recapture_std < 200)) +
  			         geom_sf(data = Norway) + theme_bw() +
  			  geom_jitter(aes(x=cLon, y=cLat, col=ID),size=1, width=0.1, height=0.1) +
  			  geom_point(aes(x=Longitude, y=Latitude, col=ID),size=1) +
  			  scale_color_viridis_c()


        datplot <- Data_mackerel_use_Ireland_select %>% filter(julian_recapture_std < 230)
        datplot <- datplot %>% mutate(Recapture_timing = ifelse(julian_recapture_std < 200, "early",
                                                                ifelse(julian_recapture_std >= 200 &julian_recapture_std < 210, "mid1",
                                                                       ifelse(julian_recapture_std >= 210 &julian_recapture_std < 220, "mid2", "late"))))
        datplot$Recapture_timing <- as.factor(datplot$Recapture_timing)
        datplot$Recapture_timing <- factor(datplot$Recapture_timing, levels=c("early","mid1","mid2","late"))
        datplot$Catch_year <- as.numeric(as.character(datplot$Catch_year))

  			# Area effect on early recapture date
  			table1 <- datplot %>% group_by(Tag_area, Recapture_timing) %>%
  			  summarize(n=n(), duration=round(mean(duration),1), release = round(mean(julian_release_std),1))
  			table1$Longitude = c(-18,-18,-18,-18,-20,-20,-20,-20)
  			table1$Latitude= c(56.5,56.5,56.5,56.5,53,53,53,53)
  			table1$label= apply(table1, 1, function(x) paste0("n=",x[3],"\n(",x[4],"days)"))

  			label_names <-as_labeller(c(
  			  'early'="Recapture <200 julian days",
  			  'mid1'="Recapture [200;210) julian days",
  			  'mid2'="Recapture [210;220) julian days",
  			  'late'="Recapture [220;230) julian days"
  			))

  			g1 <- ggplot(data=datplot) +
  			         geom_sf(data = Norway) + theme_bw() +
  			  geom_jitter(aes(x=cLon, y=cLat, col=Catch_year, size=Catch_year), width=0.2, height=0.2, alpha=0.5) +
  			  geom_jitter(aes(x=Longitude, y=Latitude, col=Catch_year, size=Catch_year), width=0.2, height=0.2, alpha=0.5) +
  			  geom_label(data=table1, aes(x=Longitude, y=Latitude, label=label)) +
  			  # geom_text(data=title, aes(x=x,y=y,label=label)) +
  			  facet_wrap(~Recapture_timing, ncol=2, labeller= label_names) +
  			  scale_color_viridis_c(name="Year", breaks=seq(2014, 2020, by=1)) +
  			  scale_size_continuous(name = "Year", breaks=seq(2014, 2020, by=1), range=c(1,4)) +
  			  labs(x="Longitude", y="Latitude")+
  			  guides(color= guide_legend(), size=guide_legend()) +
  			theme(axis.title = element_text(size=14),
  			      axis.text = element_text(size=13),
  			      text = element_text(size=12),
  			      legend.title = element_text(size=15, hjust=0.5),
  			      legend.text = element_text(size=13),
  			      strip.text = element_text(size = 14),
  			      plot.title = element_text(hjust = 0.5, size=17))

  			  ggsave(g1, filename = "MS/figs/Recap/Iceland_recap.pdf",
  			       width=10, height=8, units="in", dpi = 450)



  			# Release timing effect on early recapture date
  			table2 <- datplot %>% group_by(Release_timing, Recapture_timing) %>%
  			  summarize(n=n(), duration=round(mean(duration),1), release = round(mean(julian_release_std),1),
  			            mvt_rate = mean(log_rate))
  			table2$Longitude = c(-24,-24,-24,-24,-24,-24,-24,-24)
  			table2$Latitude= c(56.5,56.5,56.5,56.5,53,53,53,53)
  			table2$label= apply(table2, 1, function(x) paste0(x[1]," n=",x[3],"\n(",x[4],"days)"))

  			label_names <-as_labeller(c(
  			  'early'="Recapture <200 julian days",
  			  'mid1'="Recapture [200;210) julian days",
  			  'mid2'="Recapture [210;220) julian days",
  			  'late'="Recapture [220;230) julian days"
  			))

  			g2 <- ggplot(data=datplot) +
  			         geom_sf(data = Norway) + theme_bw() +
  			  geom_jitter(aes(x=cLon, y=cLat, col=Release_timing, size=Catch_year), width=0.2, height=0.2, alpha=0.5) +
  			  geom_jitter(aes(x=Longitude, y=Latitude, col=Release_timing, size=Catch_year), width=0.2, height=0.2, alpha=0.5) +
  			  geom_label(data=table2, aes(x=Longitude, y=Latitude, label=label)) +
  			  # geom_text(data=title, aes(x=x,y=y,label=label)) +
  			  facet_wrap(~Recapture_timing, ncol=2, labeller= label_names) +
  			  # scale_color_viridis_c(name="Year", breaks=seq(2014, 2020, by=1)) +
  			  scale_size_continuous(name = "Year", breaks=seq(2014, 2020, by=1), range=c(1,4)) +
  			  labs(x="Longitude", y="Latitude")+
  			  guides(color= guide_legend(), size=guide_legend()) +
  			theme(axis.title = element_text(size=14),
  			      axis.text = element_text(size=13),
  			      text = element_text(size=12),
  			      legend.title = element_text(size=15, hjust=0.5),
  			      legend.text = element_text(size=13),
  			      strip.text = element_text(size = 14),
  			      plot.title = element_text(hjust = 0.5, size=17))

  			  ggsave(g2, filename = "MS/figs/Recap/Iceland_recap_releasetiming.pdf",
  			       width=10, height=8, units="in", dpi = 450)


  			  # Mackerel size effect on early recapture date
  			  table3 <- datplot %>% group_by(length_bin, Recapture_timing) %>%
  			    summarize(n=n(), duration=round(mean(duration),1), release = round(mean(julian_release_std),1),
  			              mvt_rate = mean(log_rate))
  			  table3$Longitude = c(-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24)
  			  table3$Latitude= c(59,59,59,59,56,56,56,56,53,53,53,53)
  			  table3$label= apply(table3, 1, function(x) paste0(x[1]," n=",x[3],"\n(",x[4],"days)"))

  			  label_names <-as_labeller(c(
  			    'early'="Recapture <200 julian days",
  			    'mid1'="Recapture [200;210) julian days",
  			    'mid2'="Recapture [210;220) julian days",
  			    'late'="Recapture [220;230) julian days"
  			  ))

  			  g3 <- ggplot(data=datplot) +
  			    geom_sf(data = Norway) + theme_bw() +
  			    geom_jitter(aes(x=cLon, y=cLat, col=length_bin, size=Catch_year), width=0.2, height=0.2, alpha=0.5) +
  			    geom_jitter(aes(x=Longitude, y=Latitude, col=length_bin, size=Catch_year), width=0.2, height=0.2, alpha=0.5) +
  			    geom_label(data=table3, aes(x=Longitude, y=Latitude, label=label)) +
  			    # geom_text(data=title, aes(x=x,y=y,label=label)) +
  			    facet_wrap(~Recapture_timing, ncol=2, labeller= label_names) +
  			    # scale_color_viridis_c(name="Year", breaks=seq(2014, 2020, by=1)) +
  			    scale_size_continuous(name = "Year", breaks=seq(2014, 2020, by=1), range=c(1,4)) +
  			    labs(x="Longitude", y="Latitude")+
  			    guides(color= guide_legend(), size=guide_legend()) +
  			    theme(axis.title = element_text(size=14),
  			          axis.text = element_text(size=13),
  			          text = element_text(size=12),
  			          legend.title = element_text(size=15, hjust=0.5),
  			          legend.text = element_text(size=13),
  			          strip.text = element_text(size = 14),
  			          plot.title = element_text(hjust = 0.5, size=17))

  			  ggsave(g3, filename = "MS/figs/Recap/Iceland_recap_fishsize.pdf",
  			         width=10, height=8, units="in", dpi = 450)


  			ggmap(map_area) + geom_jitter(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat, col=duration, pch=Tag_area)) +
  			  scale_color_viridis_c()

  			ggmap(map_area) + geom_jitter(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat, col=Latitude, cex=log_rate)) +
  			  scale_color_viridis_c()

  			with(subset(Data_mackerel_use_Ireland_select, cLon< -9), table(Tag_area, Release_year))


        # no threshold model
			    if (model_choice == "lm"){
			      use_version_simple <- paste0(getwd(), "/src/mackerel_mvt_model_nothresh")
            compile(paste0(use_version_simple, ".cpp"))
            dyn.load(use_version_simple)

            m1 <- lm(log_rate ~ factor(Tag_area)*factor(Release_timing) + factor(Release_year) + length_scaled + julian_recapture_scaled, data=Data_mackerel_use_Ireland_select)
            m_frame <- model.frame(m1)
            XX <- model.matrix(m1, m_frame)

            N_threshold <- 1
            data_tmb_nothres <- list(N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
                             X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
                             Nthres=length(threshold_vals),
                             y = Data_mackerel_use_Ireland_select$log_rate,
                             X_pred =as.matrix(as.data.frame(XX)),
                             N_pred =nrow(Data_mackerel_use_Ireland_select),
                             Likconfig = 0      # 0 = dnorm, 1 = dgamma
            )

            parameters_tmb_nothres <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                                   log_sigma = rep(log(0.2),N_threshold)
            )

            Map = list()
            obj <- MakeADFun(data_tmb_nothres, parameters_tmb_nothres, random = NULL, DLL = "mackerel_mvt_model_nothresh", map=Map)
            opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
            opt$objective * 2 + 2 * length(opt$par)

            AIC_nothres <- c()
            Maps_nothres <- matrix(1:(ncol(XX)), ncol = (ncol(XX)), nrow =ncol(XX)+1,byrow=TRUE)
            for(k in 1:(nrow(Maps_nothres)-1)){
                Map = list()
                Map$beta <- factor(Maps_nothres[k,])
                Map$log_sigma <- factor(c(ncol(XX)+1))

                #-- optimize-- :
                obj <- MakeADFun(data_tmb_nothres, parameters_tmb_nothres, random = NULL, DLL = "mackerel_mvt_model_nothresh", map=Map)
                opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                                      control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

                map <- as.numeric(c(as.character(Map$beta), (length(opt$par))))
                AIC_nothres[k] <- opt$objective * 2 + 2 * length(opt$par)
                # only removing 1
                mle <- sapply(1:length(Map$beta), function(x) ifelse(is.na(x)==F, opt$par[Map$beta[x]], NA))
                cov <- solve(obj$he(opt$par))
                sde <- sqrt(diag(cov))
                pvalue_diff <- function(i){
                  if (! TRUE %in% mle[i,]) out <- 2*pnorm(abs(diff(mle[i,])), 0, sqrt(cov[i,i]+cov[i+ncol(XX),i+ncol(XX)]-2*cov[i,i+ncol(XX)]), lower.tail = FALSE)
                  if (TRUE %in% mle[i,]) out <- NA
                  return(out)
                }
                j <- which.max(2*pnorm(abs(opt$par)/sde, lower.tail = FALSE)[-(length(opt$par))])
                j.corrected <- as.numeric(map[!is.na(map) & !duplicated(map)][j])
                Maps_nothres[(k+1):nrow(Maps_nothres), j.corrected] <- NA_real_
                parameters_tmb_nothres$beta[ifelse(j<=ncol(XX),j,j-ncol(XX)),ifelse(j/ncol(XX)>1,2,1)] <- 0
              }

              best_mod <- which.min(AIC_nothres)
              Map = list()
              Map$beta <- factor(Maps_nothres[best_mod,])
              Map$log_sigma <- factor(c(2*ncol(XX)+1))
              parameters_tmb_nothres <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                                     log_sigma = rep(log(0.2),N_threshold))
              parameters_tmb_nothres$beta[which(is.na(Maps_nothres[best_mod,]))] <- 0

              obj <- MakeADFun(data_tmb_nothres, parameters_tmb_nothres, random = NULL, DLL = "mackerel_mvt_model_nothresh", map=Map)
              opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                              control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))


              sd_report <- sdreport(obj)
              check_estimability(obj)
              sigma <- as.vector(exp(summary(sd_report, "fixed")[grep("log_sigma", rownames(summary(sd_report, "fixed"))),1]))
              mu_pred <- matrix(summary(sd_report, "report")[grep("mu_pred", rownames(summary(sd_report, "report"))),1], ncol=1, byrow=FALSE)

              par(mfrow=c(2,1), mar=c(4,3,1,1), oma=c(1,1,1,1))
              plot(mu_pred, data_tmb_nothres$y); abline(0,1)
              qqnorm(y=(mu_pred-data_tmb_nothres$y)/sd(mu_pred-data_tmb_nothres$y), xlim=c(-3.5,3.5),ylim=c(-3.5,3.5))
              abline(0,1, lty=2)

               ## Now producing the residual vs predictor
              opt_resid <- (mu_pred-data_tmb_nothres$y)/sd(mu_pred-data_tmb_nothres$y)
              par(mfrow=c(3,2))
              plot(Data_mackerel_use_Ireland_select$Length, opt_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Tag_area), opt_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Release_timing), opt_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$Release_year, opt_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$julian_recapture_std, opt_resid); abline(h=0, lty=2)

			    }

			    if (model_choice == "lme"){
			      use_version_simple <- paste0(getwd(), "/src/mackerel_mvt_model_nothresh_RE")
            compile(paste0(use_version_simple, ".cpp"))
            dyn.load(use_version_simple)

            m1 <- lm(log_rate ~ factor(Tag_area)*factor(Release_timing) + length_scaled + julian_recapture_scaled, data=Data_mackerel_use_Ireland_select)
            m1 <- lm(log_rate ~ factor(Release_timing) + Latitude_scaled:factor(Release_timing) + Latitude2_scaled:factor(Release_timing) + length_scaled + julian_recapture_scaled, data=Data_mackerel_use_Ireland_select)
            m_frame <- model.frame(m1)
            XX <- model.matrix(m1, m_frame)

            N_threshold <- 1
            data_tmb_nothres <- list(N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
                             X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
                             Nthres=length(threshold_vals),
                             y = Data_mackerel_use_Ireland_select$log_rate,
                             X_pred =as.matrix(as.data.frame(XX)),
                             N_pred =nrow(Data_mackerel_use_Ireland_select),
                             N_year = length(unique(Data_mackerel_use_Ireland_select$Release_year)),
                             Year_ID = as.numeric(as.character(Data_mackerel_use_Ireland_select$Release_year))-2014,
                             Year_ID_pred = as.numeric(as.character(Data_mackerel_use_Ireland_select$Release_year))-2014,
                             Likconfig = 0      # 0 = dnorm, 1 = dgamma
            )

            parameters_tmb_nothres <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                                   log_sigma = log(0.2),
                                   year = rep(0, length(unique(Data_mackerel_use_Ireland_select$Release_year))),
                                   log_sigma_year = log(0.2)
            )

            Map = list()

            obj <- MakeADFun(data_tmb_nothres, parameters_tmb_nothres, random = "year", DLL = "mackerel_mvt_model_nothresh_RE", map=Map)
            opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
            opt$objective * 2 + 2 * length(opt$par)

            AIC_nothres <- c()
            Maps_nothres <- matrix(1:(ncol(XX)), ncol = (ncol(XX)), nrow =ncol(XX)+1,byrow=TRUE)
            for(k in 1:(nrow(Maps_nothres)-1)){
                Map = list()
                Map$beta <- factor(Maps_nothres[k,])
                Map$log_sigma <- factor(c(ncol(XX)+1))
                Map$year <- factor(ncol(XX)+1+1:length(parameters_tmb_nothres$year))
                Map$log_sigma_year <- factor(ncol(XX)+8)

                #-- optimize-- :
                obj <- MakeADFun(data_tmb_nothres, parameters_tmb_nothres, random = "year", DLL = "mackerel_mvt_model_nothresh_RE", map=Map)
                opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                                      control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

                map <- as.numeric(c(as.character(Map$beta), (length(opt$par)-c(1,0))))
                AIC_nothres[k] <- opt$objective * 2 + 2 * length(opt$par)
                # only removing 1
                mle <- sapply(1:length(Map$beta), function(x) ifelse(is.na(x)==F, opt$par[Map$beta[x]], NA))
                sdrep <- sdreport(obj)
                cov <- sdrep$cov.fixed
                sde <- summary(sdrep, "fixed")[,2]
                pvalue_diff <- function(i){
                  if (! TRUE %in% mle[i,]) out <- 2*pnorm(abs(diff(mle[i,])), 0, sqrt(cov[i,i]+cov[i+ncol(XX),i+ncol(XX)]-2*cov[i,i+ncol(XX)]), lower.tail = FALSE)
                  if (TRUE %in% mle[i,]) out <- NA
                  return(out)
                }
                j <- which.max(2*pnorm(abs(opt$par)/sde, lower.tail = FALSE)[-(length(opt$par)-0:1)])
                j.corrected <- as.numeric(map[!is.na(map) & !duplicated(map)][j])
                Maps_nothres[(k+1):nrow(Maps_nothres), j.corrected] <- NA_real_
                parameters_tmb_nothres$beta[ifelse(j<=ncol(XX),j,j-ncol(XX)),ifelse(j/ncol(XX)>1,2,1)] <- 0
              }

              best_mod <- which.min(AIC_nothres)
              Map = list()
              Map$beta <- factor(Maps_nothres[best_mod,])
              Map$log_sigma <- factor(c(2*ncol(XX)+1))
              Map$year <- factor(ncol(XX)+1+1:length(parameters_tmb_nothres$year))
              Map$log_sigma_year <- factor(ncol(XX)+8)
              parameters_tmb_nothres <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                                     log_sigma = log(0.2),
                                     year = rep(0, length(unique(Data_mackerel_use_Ireland_select$Release_year))),
                                     log_sigma_year = log(0.2)
              )
              parameters_tmb_nothres$beta[which(is.na(Maps_nothres[best_mod,]))] <- 0

              obj <- MakeADFun(data_tmb_nothres, parameters_tmb_nothres, random = "year", DLL = "mackerel_mvt_model_nothresh_RE", map=Map)
              opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                              control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

              AIC_nothres_best <- opt$objective * 2 + 2 * length(opt$par)
              sd_report_nothres <- sdreport(obj)
              check_estimability(obj)
              sigma_nothres <- as.vector(exp(summary(sd_report_nothres, "fixed")[grep("log_sigma", rownames(summary(sd_report_nothres, "fixed"))),1]))
              mu_pred_nothres <- matrix(summary(sd_report_nothres, "report")[grep("mu_pred", rownames(summary(sd_report_nothres, "report"))),1], ncol=1, byrow=FALSE)

              par(mfrow=c(2,1), mar=c(4,3,1,1), oma=c(1,1,1,1))
              plot(mu_pred_nothres, data_tmb_nothres$y); abline(0,1)
              qqnorm(y=(mu_pred_nothres-data_tmb_nothres$y)/sd(mu_pred_nothres-data_tmb_nothres$y))
              abline(0,1, lty=2)

               ## Now producing the residual vs predictor
              opt_resid <- (mu_pred_nothres-data_tmb_nothres$y)/sd(mu_pred_nothres-data_tmb_nothres$y)
              d <- Data_mackerel_use_Ireland_select
              d$residuals <- opt_resid
              d <- d[order(d$residuals),]
              d$qq <- (qnorm(ppoints(nrow(d))))
              p1 <- ggplot(d, aes(x= qq, y=residuals, col=julian_recapture_std)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw() +
                scale_color_viridis_c(option = "plasma") + labs(x="Theoretical", y="Sample")
              p1
              par(mfrow=c(3,2))
              plot(Data_mackerel_use_Ireland_select$Length, opt_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Tag_area), opt_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Release_timing), opt_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$Release_year, opt_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$julian_recapture_std, opt_resid); abline(h=0, lty=2)

			    }


        # A single change point model
    				# model selection for the single break point model
    				  if (model_choice == "lm") source(paste0(getwd(), "/R/model_selection.R"))
    				  if (model_choice == "lme") source(paste0(getwd(), "/R/model_selection_RE.R"))

      				sd_report_1break_best <- sdreport(obj1break_best)
      				check_estimability(obj1break_best)
      				sigma <- as.vector(exp(summary(sd_report_1break_best, "fixed")[grep("log_sigma", rownames(summary(sd_report_1break_best, "fixed"))),1]))
      				mu_pred_best <- matrix(summary(sd_report_1break_best, "report")[which(rownames(summary(sd_report_1break_best, "report")) == "mu_pred"),1], ncol=2, byrow=FALSE)

    				# Calculating the actual prediction
    				  LL <- obj1break_best$report()$LL

    			  # The weighting factor is the average likelihood across the N datapoint per threshold value
    					weight <- apply(exp(LL), 2, mean)
    					weight <- weight/sum(weight)
    					plot(data_tmb$thresh, weight, type="p", pch=".",  xlab="Recapture date (Julian days)", ylab="Likelihood", cex=1.4)
    					lines(data_tmb$thresh, weight, cex=0.9, col=grey(.5))
    #           weight <- rep(0,17)
    # 					weight[5]=1

    					Prediction <- rep(0, data_tmb$N)
    					for (n in 1:data_tmb$N){
    					  for (thr in 1:data_tmb$Nthres){
    					    if (data_tmb$thres_cov[n] < (data_tmb$thresh[thr])){
    					      Prediction[n] = Prediction[n] + weight[thr]*mu_pred_best[n,1]
    					    }
    					    if (data_tmb$thres_cov[n] >= (data_tmb$thresh[thr])){
    					      Prediction[n] = Prediction[n] + weight[thr]*mu_pred_best[n,2]
    					    }
    					  }
    					}


    					# Calculating the R2 (nakagawa's approach)
    					mu_pred_fixed <- matrix(summary(sd_report_1break_best, "report")[grep("mu_pred_fixed", rownames(summary(sd_report_1break_best, "report"))),1], ncol=2, byrow=FALSE)
    					Prediction_fixed <- rep(0, data_tmb$N)
    					for (n in 1:data_tmb$N){
    					  for (thr in 1:data_tmb$Nthres){
    					    if (data_tmb$thres_cov[n] < (data_tmb$thresh[thr])){
    					      Prediction_fixed[n] = Prediction_fixed[n] + weight[thr]*mu_pred_fixed[n,1]
    					    }
    					    if (data_tmb$thres_cov[n] >= (data_tmb$thresh[thr])){
    					      Prediction_fixed[n] = Prediction_fixed[n] + weight[thr]*mu_pred_fixed[n,2]
    					    }
    					  }
    					}
    					var_fixed <- var(Prediction_fixed-data_tmb$y)
    					var_obs <- obj1break_best$report()$var_obs
    					var_year <- obj1break_best$report()$var_year

    					var(mu_pred_fixed[,1])/(var_obs[1]+var_year[1]+var(mu_pred_fixed[,1]))
    					var(mu_pred_fixed[,2])/(var_obs[2]+var_year[2]+var(mu_pred_fixed[,2]))


    					## Now producing the residual vs predictor

    					opt1break_resid <- (Prediction-data_tmb$y)/sd(Prediction-data_tmb$y)
              qqnorm(opt1break_resid); abline(0,1)
    					par(mfrow=c(3,2))
              plot(Data_mackerel_use_Ireland_select$Length, opt1break_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Tag_area), opt1break_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Release_timing), opt1break_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$Release_year, opt1break_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$julian_recapture_std, opt1break_resid); abline(h=0, lty=2)


              d <- Data_mackerel_use_Ireland_select
              d$residuals <- opt1break_resid
              d <- d[order(d$residuals),]
              d$qq <- (qnorm(ppoints(nrow(d))))
              p1 <- ggplot(d, aes(x= qq, y=residuals, col=julian_recapture_std)) + geom_point(size=1) + geom_abline(slope=1) + theme_bw() +
                scale_color_viridis_c(option = "plasma") + labs(x="Theoretical", y="Sample")
              p1

              ## main effect table & plot
              map <- as.factor(c(as.numeric(as.character(Map_best$beta)), max(as.numeric(as.character(Map_best$beta)),na.rm=T)+1:length(grep("sigma", names(opt1break_best$par)))))
              map <- as.numeric(factor(map, labels=1:length(levels(map))))
              mle <- sapply(1:length(map), function(x) ifelse(is.na(x)==F, opt1break_best$par[map[x]], NA))
              cov <- sd_report_1break_best$cov.fixed
              sde <- sapply(1:length(map), function(x) ifelse(is.na(x)==F, sqrt(diag(cov))[map[x]], NA))
              Estimates <- cbind(mle, sde)
              Estimates <- rbind(Estimates, summary(sd_report_1break_best, "random"))
              rownames(Estimates) <- c(apply(cbind(colnames(XX),"before"),1, function(x) paste(x, collapse="_")), apply(cbind(colnames(XX),"after"),1, function(x) paste(x, collapse="_")),
                                       apply(cbind(names(opt1break_best$par)[grep("sigma", names(opt1break_best$par))], rep(c("before", "after"),2)), 1, function(x) paste(x,collapse="_")),
                                       apply(expand.grid(year=2014:2020, time=c("before", "after")), 1, function(x) paste(x,collapse="_")))
              Estimates <- as.data.frame(Estimates)
              #Estimates$gradient <- rep(NA, length(map))
              #Estimates$gradient <- opt1break_best$diagnostics$final_gradient[map[!is.na(map)]]

              fake <- data.frame(Threshold = data_tmb$thresh, Likelihood = weight)
              plot(data_tmb$thresh, weight, type="p", main="", xlab="", ylab="", pch=".", cex=1.5)
              lines(data_tmb$thresh, weight, cex=0.9, col=grey(0.5))
              mtext(side = 1, line=3, "Recapture date (julian days)", cex=1.2)
              mtext(side = 2, line=3, "Likelihood", cex=1.2)
              #mtext(side = 3, line=1, "Threshold value", lwd=2, cex=1.4)
              abline(v=243, lty=2)
              abline(v=273, lty=2)
              abline(v=304, lty=2)
              # text(x=225, y=0.057, "August")
              text(x=258, y=0.057, "September")
              text(x=288, y=0.057, "October")

              g_threshold <- ggplot(data=fake, aes(x=Threshold, y=Likelihood)) + geom_line() + theme_bw() + geom_point(pch=1, size=0.3)  +
                geom_vline(xintercept=240, linetype=2) + geom_vline(xintercept=270, linetype=2) + geom_vline(xintercept=300, linetype=2) +
                geom_text(data=data.frame(Threshold=255, Likelihood=max(fake$Likelihood)*1.02), label="September", size=5) +
                geom_text(data=data.frame(Threshold=285, Likelihood=max(fake$Likelihood)*1.02), label="October", size=5) +
                #labs(x="Recapture time (julian days)")
                labs(x="Recapture time (julian days)", title="Threshold value") +
                theme(axis.title = element_text(size=15)) + theme(axis.title = element_text(size=15),
                                                                  axis.text = element_text(size=13),
                                                                  plot.title = element_text(hjust = 0.5, size=20))
              g_threshold
              ggsave(g_threshold, filename = "MS/figs/Fig5.pdf",
                     width=18, height=13, units="cm", dpi = 450)


              Estimates_df <- as.data.frame(Estimates[,1:2])

              colnames(Estimates_df) <- c("mean", "se")
              Estimates_df$group <- as.character(sapply(as.character(sapply(rownames(Estimates_df), function(x) str_sub(x, start=-6))), function(x) str_remove(x, "_")))
              Estimates_df$group <- factor(Estimates_df$group, levels=c("before","after"))
              Estimates_df$group <- factor(Estimates_df$group, labels=c("Before threshold","After threshold"))
              # Estimates_df$covariate <- c(rep(c("Intercept", "From_South", "> May 22nd", "body length (rescaled)", "recapture date (rescaled)", "From South + >May 22nd"), 2),
              #                             "log(sigma_obs)", "log(sigma_obs)", "log(sigma_year)", "log(sigma_year)", rep(2014:2020, 2))
              Estimates_df$covariate <- c(rep(c("Intercept", ">= May 22nd", "body length (rescaled)", "recapture date (rescaled)", "Latitude + < May 22nd", "Latitude + >= May 22nd"
                                                , "Latitude^2 + < May 22nd", "Latitude^2 + >= May 22nd"), 2),
                                          "log(sigma_obs)", "log(sigma_obs)", "log(sigma_year)", "log(sigma_year)", rep(2014:2020, 2))
              # Estimates_df$covariate <- c(rep(c("Intercept", "> May 22nd", "body length (rescaled)", "recapture date (rescaled)", "Latitude", "Latitude^2"),2),
              #                              "log(sigma_obs)", "log(sigma_obs)", "log(sigma_year)", "log(sigma_year)", rep(2014:2020, 2))
              # #Estimates_df$mean[grep("sigma", Estimates_df$covariate)] <- exp(Estimates_df[grep("sigma", Estimates_df$covariate),'mean'])
              #Estimates_df$covariate <- c(rep(c("Intercept", "From_South", "> May 22nd", "body length (rescaled)", "recapture date (rescaled)", "From South + >May 22nd"), 2),
              #                            "sigma_obs", "sigma_obs", "sigma_year", "sigma_year")
              #Estimates_df$covariate <- factor(Estimates_df$covariate, levels=c("Intercept", "> May 22nd", "From South + >May 22nd", "body length (rescaled)", "recapture date (rescaled)", "log(sigma_obs)", "log(sigma_year)", 2014:2020))
              # Estimates_df$covariate <- factor(Estimates_df$covariate, levels=c("Intercept", "> May 22nd", "Latitude", "Latitude^2", "body length (rescaled)", "recapture date (rescaled)", "log(sigma_obs)", "log(sigma_year)", 2014:2020))
              Estimates_df$covariate <- factor(Estimates_df$covariate, levels=c("Intercept", ">= May 22nd", "Latitude + < May 22nd", "Latitude + >= May 22nd", "Latitude^2 + < May 22nd", "Latitude^2 + >= May 22nd", "body length (rescaled)", "recapture date (rescaled)", "log(sigma_obs)", "log(sigma_year)", 2014:2020))
              levels(Estimates_df$covariate)[1] <- "Intercept (true_value add 8.0)"
              Estimates_df$mean_transformed <- ifelse(Estimates_df$mean>8, Estimates_df$mean-8, Estimates_df$mean)
              Estimates_df$x_fake = 3.5
              Estimates_df$Label = paste0("", round(Estimates_df$mean, 2), "(", round(Estimates_df$se, 3), ")")
              Estimates_df <- Estimates_df[1:20,]

              library(ggplot2)
              fig6 <- ggplot(Estimates_df, aes(x=mean_transformed, y=covariate)) + facet_grid(.~ group) + geom_point() +
                geom_errorbar(aes(xmin=mean_transformed-1.96*se, xmax=mean_transformed+1.96*se, width=0.1))+
                geom_vline(xintercept=0, linetype=2) + theme_bw() + coord_cartesian(xlim=c(-4,6.3)) +
                geom_text(aes(x=x_fake, y=covariate, label=Label),size=3,hjust=0) +
                labs(x="Values", y="Covariate") +
                theme(axis.title = element_text(size=15),
                      axis.text = element_text(size=12),
                      legend.text = element_text(size=12),
                      strip.text = element_text(size = 15),
                      plot.title = element_text(hjust = 0.5, size=20))
              ggsave(fig6, filename = "MS/figs/Fig6.pdf",
                     width=8.27, height=5.83, units="in", dpi = 450)


              # Marginal effect of release latitude
              latitude <- seq(min(Data_mackerel_use_Ireland_select$Latitude), max(Data_mackerel_use_Ireland_select$Latitude),by=0.01)
              latitude_scaled <- (latitude-mean(Data_mackerel_use_Ireland_select$Latitude))/sd(Data_mackerel_use_Ireland_select$Latitude)
              latitude2 <- latitude^2
              latitude2_scaled <- (latitude2 - mean((Data_mackerel_use_Ireland_select$Latitude)^2))/sd((Data_mackerel_use_Ireland_select$Latitude)^2)

              par(mfrow=c(2,2), mar=c(4,2,2,2), oma=c(1,3,0,0), cex.axis=1.2, cex.lab=1.3, cex.main=1.4)
              # the effect of latitude (release location)
                lat_effect_before <- latitude_scaled*Estimates_df[5,1]+Estimates_df[7,1]*latitude2_scaled
                lat_effect_before_low <- latitude_scaled*Estimates_df[5,1]+Estimates_df[7,1]*latitude2_scaled - sqrt(latitude_scaled^2*Estimates_df[5,2]^2+Estimates_df[7,2]^2*latitude2_scaled^2)
                lat_effect_before_high <- latitude_scaled*Estimates_df[5,1]+Estimates_df[7,1]*latitude2_scaled + sqrt(latitude_scaled^2*Estimates_df[5,2]^2+Estimates_df[7,2]^2*latitude2_scaled^2)
                plot(latitude, lat_effect_before, type="l", lwd=2, xlab="Latitude (°)", ylab="Marginal effect", main ="", ylim=c(-0.8,0.55))
                polygon(c(latitude, rev(latitude)), c(lat_effect_before_low, rev(lat_effect_before_high)), col=grey(0.5,alpha=0.5), border=grey(0.5,alpha=0.5))
                lat_effect_after <- latitude_scaled*Estimates_df[13,1]+Estimates_df[15,1]*latitude2_scaled
                lat_effect_after_low <- latitude_scaled*Estimates_df[13,1]+Estimates_df[15,1]*latitude2_scaled - sqrt(latitude_scaled^2*Estimates_df[13,2]^2+Estimates_df[15,2]^2*latitude2_scaled^2)
                lat_effect_after_high <- latitude_scaled*Estimates_df[13,1]+Estimates_df[15,1]*latitude2_scaled + sqrt(latitude_scaled^2*Estimates_df[13,2]^2+Estimates_df[15,2]^2*latitude2_scaled^2)
                polygon(c(latitude, rev(latitude)), c(lat_effect_after_low, rev(lat_effect_after_high)), col=rgb(0.8,0,0,alpha=0.3), border=rgb(0.8,0,0,alpha=0.3))
                rug(Data_mackerel_use_Ireland_select$Latitude)

                lat_effect_before2 <- latitude_scaled*Estimates_df[6,1]+Estimates_df[8,1]*latitude2_scaled
                lat_effect_before2_low <- latitude_scaled*Estimates_df[6,1]+Estimates_df[8,1]*latitude2_scaled - sqrt(latitude_scaled^2*Estimates_df[6,2]^2+Estimates_df[8,2]^2*latitude2_scaled^2)
                lat_effect_before2_high <- latitude_scaled*Estimates_df[6,1]+Estimates_df[8,1]*latitude2_scaled + sqrt(latitude_scaled^2*Estimates_df[6,2]^2+Estimates_df[8,2]^2*latitude2_scaled^2)
                polygon(c(latitude, rev(latitude)), c(lat_effect_before2_low, rev(lat_effect_before2_high)), col=grey(0.2,alpha=0.5), border=grey(0.2,alpha=0.5))
                lat_effect_after2 <- latitude_scaled*Estimates_df[14,1]+Estimates_df[16,1]*latitude2_scaled
                lat_effect_after2_low <- latitude_scaled*Estimates_df[14,1]+Estimates_df[16,1]*latitude2_scaled - sqrt(latitude_scaled^2*Estimates_df[14,2]^2+Estimates_df[16,2]^2*latitude2_scaled^2)
                lat_effect_after2_high <- latitude_scaled*Estimates_df[14,1]+Estimates_df[16,1]*latitude2_scaled + sqrt(latitude_scaled^2*Estimates_df[14,2]^2+Estimates_df[16,2]^2*latitude2_scaled^2)
                polygon(c(latitude, rev(latitude)), c(lat_effect_after2_low, rev(lat_effect_after2_high)), col=rgb(1,0,0,alpha=0.3), border=rgb(1,0,0,alpha=0.3))

                lines(latitude, lat_effect_before, col="black", lwd=2)
                lines(latitude, lat_effect_after, col="red", lwd=2)
                lines(latitude, lat_effect_before2, col="black", lwd=2, lty=2)
                lines(latitude, lat_effect_after2, col="red", lwd=2, lty=2)
                legend("bottomleft", c("<May 22nd", expression("">=May22nd)), col=c("black", "black"), lty=c(1,2), lwd=2, bty="n")

              # the effect of recapture date
                reca <- seq(min(Data_mackerel_use_Ireland_select$julian_recapture_std), max(Data_mackerel_use_Ireland_select$julian_recapture_std), by=0.01)
                reca_scaled <- (reca-mean(Data_mackerel_use_Ireland_select$julian_recapture_std))/sd(Data_mackerel_use_Ireland_select$julian_recapture_std)
                recap_effect_before <- reca_scaled*Estimates_df[4,1]
                recap_effect_before_low <- reca_scaled*Estimates_df[4,1] - sqrt(reca_scaled^2*Estimates_df[4,2]^2)
                recap_effect_before_high <- reca_scaled*Estimates_df[4,1] + sqrt(reca_scaled^2*Estimates_df[4,2]^2)
                plot(reca, recap_effect_before, type="l", lwd=2, xlab="Recapture date (julian days)", ylab="Marginal effect", main ="", ylim=c(-1, 1.80))
                polygon(c(reca, rev(reca)), c(recap_effect_before_low, rev(recap_effect_before_high)), col=grey(0.5, alpha=0.5), border=grey(0.5, alpha=0.5))
                recap_effect_after <- reca_scaled*Estimates_df[12,1]
                recap_effect_after_low <- reca_scaled*Estimates_df[12,1] - sqrt(reca_scaled^2*Estimates_df[12,2]^2)
                recap_effect_after_high <- reca_scaled*Estimates_df[12,1]+ sqrt(reca_scaled^2*Estimates_df[12,2]^2)
                polygon(c(reca, rev(reca)), c(recap_effect_after_low, rev(recap_effect_after_high)), col=rgb(1,0,0,alpha=0.3), border=rgb(1,0,0,alpha=0.3))
                lines(reca, recap_effect_before, col="black", lwd=2)
                lines(reca, recap_effect_after, col="red", lwd=2)
                rug(Data_mackerel_use_Ireland_select$julian_recapture_std)
                # abline(v=212, lty=2)
                # abline(v=243, lty=2)
                # abline(v=273, lty=2)
                # abline(v=304, lty=2)
                # text(x=198, y=1.6, "July", srt=45)
                # text(x=228, y=1.6, "August", srt=45)
                # text(x=258, y=1.6, "September", srt=45)
                # text(x=288, y=1.6, "October", srt=45)
                # text(x=318, y=1.6, "November", srt=45)
                # the effect of fish size
                ll <- seq(min(Data_mackerel_use_Ireland_select$Length), max(Data_mackerel_use_Ireland_select$Length), by=0.01)
                ll_scaled <- (ll-mean(Data_mackerel_use_Ireland_select$Length))/sd(Data_mackerel_use_Ireland_select$Length)
                length_effect_before <- ll_scaled*Estimates_df[3,1]
                length_effect_before_low <- ll_scaled*Estimates_df[3,1] - sqrt(ll_scaled^2*Estimates_df[3,2]^2)
                length_effect_before_high <- ll_scaled*Estimates_df[3,1] + sqrt(ll_scaled^2*Estimates_df[3,2]^2)
                plot(ll, length_effect_before, type="l", lwd=2, xlab="Length (cm)", ylab="Marginal effect", main ="", ylim=c(-0.08,0.08))
                polygon(c(ll, rev(ll)), c(length_effect_before_low, rev(length_effect_before_high)), col=grey(0.5, alpha=0.5), border=grey(0.5, alpha=0.5))
                length_effect_after <- ll_scaled*Estimates_df[11,1]
                length_effect_after_low <- ll_scaled*Estimates_df[11,1] - sqrt(ll_scaled^2*Estimates_df[11,2]^2)
                length_effect_after_high <- ll_scaled*Estimates_df[11,1] + sqrt(ll_scaled^2*Estimates_df[11,2]^2)
                polygon(c(ll, rev(ll)), c(length_effect_after_low, rev(length_effect_after_high)), col=rgb(1,0,0,alpha=0.3), border=rgb(1,0,0,alpha=0.3))
                lines(ll, length_effect_before, col="black", lwd=2)
                lines(ll, length_effect_after, col="red", lwd=2)
                rug(Data_mackerel_use_Ireland_select$Length)
              # empty plot for legend
                plot(0,0,type="n",bty="n", axes=F,xlab="", ylab="")
                legend("center", c("Before threshold", "After threshold"), col=c("black", "red"), lty=1, lwd=2, bty="n", cex=1.5)
              # ylabel
                mtext(side=2, "Marginal effect on linear movement rate", outer=T, line=1.5, cex=1.3)
                dev.copy2pdf(file=paste0(getwd(), "/MS/figs/Fig7.pdf"),
                             width=8.27, height=5.83, out.type="pdf")
                dev.off()


              # Likelihood profile & residual
              par(mfrow=c(4,4), oma=c(1,1,1,1), mar=c(2,3,3,1))
              for (i in 1:length(parameters_tmb$beta))	plot(tmbprofile(obj1break_best, name = i), main=Estimates_df$covariate[i])
              # plot(0,0,type="n", bty="n", axes=F)
              # text(x=0, y=0, "max_gradient = 0.004", xpd=NA)

              ## Now plotting the residuals
              par(mfrow=c(2,1), mar=c(4,3,1,1), oma=c(1,1,1,1))
              plot(Prediction, data_tmb$y, main="Best threshold model"); abline(0,1)
              qqnorm(y=(Prediction-data_tmb$y)/sd(Prediction-data_tmb$y), xlim=c(-3.5,3.5),ylim=c(-3.5,3.5), main="Best threshold model")
              abline(0,1, lty=2)

              ## Residuals performance comparison between the non-threshhold and threshold model
              par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,2,1,1))
              plot(mu_pred_nothres, data_tmb_nothres$y, main="Best non-threshold model", xlab="Predicted log movement rate",
                   ylab="Observed log movement rate"); abline(0,1)
              stats::qqnorm(y=(mu_pred_nothres-data_tmb_nothres$y)/sd(mu_pred_nothres-data_tmb_nothres$y), pch=20, xlim=c(-3.5,3.5),ylim=c(-3.5,3.5), main="Best non-threshold model")
              abline(0,1, lty=2)

              par(mfrow=c(1,2), mar=c(4,4,1,1), oma=c(1,2,1,1))
              plot(Prediction, data_tmb$y, main="", xlab="Predicted log movement rate",
                   ylab="Observed log movement rate"); abline(0,1)
              stats::qqnorm(y=(Prediction-data_tmb$y)/sd(Prediction-data_tmb$y), pch=20, xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),
                     xlab="Theoretical quantiles", ylab="Sample quantiles", main="")
              abline(0,1, lty=2)

              ### Do some simulations as in DHARMA to look at the model fit - DOES not work because predicted value is not ordered...
              # sims <- replicate(1000, obj1break_best$simulate()$y_sim)
              # # now making prediction using these simulated observation
              # y_sim <- matrix(0, nrow=data_tmb$N, ncol=1000)
              # for (sim in 1:1000){
              #   for (n in 1:data_tmb$N){
              #     for (thr in 1:data_tmb$Nthres){
              #       if (data_tmb$thres_cov[n] < (data_tmb$thresh[thr])){
              #         y_sim[n,sim] = y_sim[n,sim] + weight[thr]*sims[n,1,sim]
              #       }
              #       if (data_tmb$thres_cov[n] >= (data_tmb$thresh[thr])){
              #         y_sim[n,sim] = y_sim[n,sim] + weight[thr]*sims[n,2,sim]
              #       }
              #     }
              #   }
              # }
              # run_resid_check(sims_y=y_sim, yobs=data_tmb$y, pred_y=Prediction, quantiles = c(0.25,0.5,0.75), rank=TRUE)
              #
              # ijk = 100
              # hist(sims[ijk,1,])
              # hist(sims[ijk,2,])
              # hist(y_sim[ijk,]); abline(v=Prediction[ijk], lty=2, lwd=2)
              # scaled_resid <- getQuantile(y_sim, data_tmb$y, method="PIT")
              # gap::qqunif(scaled_resid,pch=2,bty="n", logscale = F, col = "black", cex = 0.6, main = "QQ uniform plot residuals", cex.main = 1)
              #


            # And now some cross-validation code to compare the two best models
              if (model_choice == "lm") source(paste0(getwd(), "/R/cross_validation.R"))
              if (model_choice == "lme") source(paste0(getwd(), "/R/cross_validation_RE.R"))


            # group together the model selection results
              library(gridExtra)
              grid.arrange(g1, g2, g_resid, g_threshold, nrow=2)


      # Now doing some prediction to illustrate the year effect along with early vs tag location

        fixed <- Estimates_df[,1:2]
        random <- summary(sd_report_1break_best, "random")
        rownames(random) <- c(paste0("Before", rep(2014:2020)), paste0("After", rep(2014:2020)))

        data_pred <- expand.grid(log_rate=1, Release_timing=c("First half", "Second half"), la_scaled=c(-0.7755247, 1.0254164), la2_scaled=c(-0.7726249, 1.0143130), length_scaled=0, julian_recapture_scaled=0)
        data_pred <- data_pred[-c(3,4,5,6),]
        mmm <- lm(log_rate ~ factor(Release_timing) + la_scaled:factor(Release_timing) + la2_scaled:factor(Release_timing) + length_scaled + julian_recapture_scaled, data=data_pred)
        mmm_frame <- model.frame(mmm)
        XXX <- model.matrix(mmm, mmm_frame)

        pred_example_before <- as.matrix(XXX)%*%fixed[1:8,1]
        bla1 <- apply((matrix(rep(pred_example_before,7), ncol=7)), 1, function(x) x+random[1:7,1])
        bla1_var <- as.matrix(XXX^2)%*%fixed[1:8,2]^2
        bla1_sd <- apply((matrix(rep(bla1_var,7), ncol=7)), 1, function(x) sqrt(x+random[1:7,2]^2))

        pred_example_after <- as.matrix(XXX)%*%fixed[9:16,1]
        bla2 <- apply((matrix(rep(pred_example_after,7), ncol=7)), 1, function(x) x+random[8:14,1])
        bla2_var <- as.matrix(XXX^2)%*%fixed[9:16,2]^2
        bla2_sd <- apply((matrix(rep(bla2_var,7), ncol=7)), 1, function(x) sqrt(x+random[8:14,2]^2))

        Before_lograte <- c()
        for (yr in seq_along(2014:2020)){
          temp1 <- data_pred
          temp1$Year <- (2014:2020)[yr]
          temp1$log_rate <- bla1[yr,]
          temp1$log_rate_sd <- bla1_sd[yr,]
          temp1$Area <- c("52N", "52N", "55N", "55N")
          temp1$Release_timing <- c("< May22nd",">= May22nd","< May22nd",">= May22nd")
          Before_lograte <- rbind(Before_lograte, temp1)
        }
        After_lograte <- c()
        for (yr in seq_along(2014:2020)){
          temp2 <- data_pred
          temp2$Year <- (2014:2020)[yr]+0.1
          temp2$log_rate <- bla2[yr,]
          temp2$log_rate_sd <- bla2_sd[yr,]
          temp2$Area <- c("52N", "52N", "55N", "55N")
          temp2$Release_timing <- c("< May22nd",">= May22nd","< May22nd",">= May22nd")
          After_lograte <- rbind(After_lograte, temp2)
        }
        Fig8 <- ggplot(Before_lograte, aes(y=log_rate, x=Year)) + geom_line(col="darkgreen") + geom_point() +
          geom_errorbar(aes(ymin=log_rate - log_rate_sd, ymax= log_rate + log_rate_sd), width=0.5, col="darkgreen") +
          geom_line(data=After_lograte, aes(y=log_rate, x=Year), col="red") + geom_point(data=After_lograte, aes(y=log_rate, x=Year)) +
          geom_errorbar(data=After_lograte, aes(ymin=log_rate - log_rate_sd, ymax= log_rate + log_rate_sd), width=0.5, col="red") +
          facet_grid(Release_timing ~ Area, labeller = labeller(type=label_parsed)) + theme_bw() +
          labs(x="Year", y="log(linear movement rate (m/day))") +
          theme(axis.title = element_text(size=12),
                axis.text = element_text(size=10),
                legend.text = element_text(size=10),
                strip.text = element_text(size = 12),
                plot.title = element_text(hjust = 0.5, size=17))
        ggsave(Fig8, filename = "MS/figs/Fig8.pdf",
               width=5.83, height=4.13, units="in", dpi = 450)






      # a bit more plot of the best threshold model
        plot(Data_mackerel_use_Ireland_select$julian_recapture_scaled, Data_mackerel_use_Ireland_select$log_rate,
             xlab="recapture data (julian days)", ylab = "log(movement rate)", cex=0.5, col=grey(0.5))
        points(subset(Data_mackerel_use_Ireland_select, Tag_area=="South_Ireland")$julian_recapture_scaled, Prediction[which(Data_mackerel_use_Ireland_select$Tag_area == "South_Ireland")], col="red3", pch=20, cex=0.5)
        points(subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland")$julian_recapture_scaled, Prediction[which(Data_mackerel_use_Ireland_select$Tag_area == "North_Ireland")], col="blue3", pch=20, cex=0.5)
        #points(Data_mackerel_use_Ireland_select$julian_recapture_scaled, mu_pred_best[,1], col="red3", pch=17, cex=0.5)
        #points(Data_mackerel_use_Ireland_select$julian_recapture_scaled, mu_pred_best[,2], col="blue3", pch=16, cex=0.5)
        dat_plot <- data.frame(Recap_date = Data_mackerel_use_Ireland_select$julian_recapture_std,
                               Before_threshold = mu_pred_best[,1],
                               After_threshold = mu_pred_best[,2],
                               Final_prediction = Prediction,
                               Obs=Data_mackerel_use_Ireland_select$log_rate)
        dat_plot_df <- reshape2::melt(dat_plot, id.vars=c("Recap_date", "Obs"), variable.name=c("Variable"))
        dat_plot_df$Variable <- as.factor(dat_plot_df$Variable)
        ggplot(dat_plot_df, aes(x=Recap_date, y=value, col=Variable, shape=Variable)) + geom_point() +
          theme_bw() + geom_point(aes(x=Recap_date, y=Obs), pch=1, col="black", cex=0.7) +
          scale_color_manual(values = c("#e69F00", "#0072B2", "#CC79A7")) #viridis_d(begin=0.1, direction = -1) +
        ggplot(subset(dat_plot_df, Variable != "Final_prediction"), aes(x=Recap_date, y=value, shape=Variable)) + geom_point(col=grey(0.99)) +
          theme_bw() + geom_point(aes(x=Recap_date, y=Obs), pch=1, col="black", cex=0.7)
        ggplot(subset(dat_plot_df, Variable != "Final_prediction"), aes(x=Recap_date, y=value, col=Variable, shape=Variable)) + geom_point() +
          theme_bw() + geom_point(aes(x=Recap_date, y=Obs), pch=1, col="black", cex=0.7) +
          scale_color_manual(values = c("#e69F00", "#0072B2")) #viridis_d(begin=0.1, direction = -1) +

        ggplot(dat_plot_df, aes(x=Recap_date, y=value, col=Variable, shape=Variable)) + geom_point() +
          theme_bw() + geom_point(aes(x=Recap_date, y=Obs), pch=16, col="black", cex=0.7) +
          scale_color_manual(values = c("#e69F00", "#0072B2", "#CC79A7")) #viridis_d(begin=0.1, direction = -1) +


        plot_movement_rate <- ggplot(Atlantic) + geom_sf() + geom_jitter(data=Data_mackerel_use_Ireland_select, aes(x=cLon, y=cLat, col=log_rate), size=2) +
          facet_wrap(~Catch_year,ncol=2) + scale_color_gradient(low="lightblue1", high="darkblue", name = " Movement \n  rate (log)") +
          labs(x="Longitude", y="Latitude") +
          theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                             legend.title = element_text(face = "bold", size = 12),
                             strip.background = element_rect(fill = "white"),
                             strip.text = element_text(face ="bold", size = 14),
                             axis.title = element_text(size=14),
                             axis.text = element_text(size=13))

        ggsave(plot_movement_rate, filename = "MS/figs/Fig_mvtrate.pdf",
               width=28, height=32, units="cm", dpi = 450)


        fig3 <- ggarrange(plot_catch, plot_movement_rate, labels=c("a)", "b)"), ncol=2)
        ggsave(fig3, filename = "MS/figs/Fig_3.pdf",
               width=36, height=20, units="cm", dpi = 450)




        # A 2 break points
            run_2breakpoint_model <- function(){
              N_threshold <- 3
              data_tmb$K <- N_threshold
              use_version_2breaks <- paste0(getwd(), "/src/mackerel_mvt_model_2breaks")
              compile(paste0(use_version_2breaks, ".cpp"))
              dyn.load(use_version_2breaks)
              #parameters_tmb <- to_save[[2]]
              set.seed(1)
              data_tmb <- list(K=N_threshold,  # number of mixture components
                               N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
                               X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
                               Nthres=length(threshold_vals),
                               thresh=threshold_vals,
                               mean_diff_tag_area= mean_diff_tag_area,
                               is_from_South=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "South_Ireland",1,0),
                               thres_cov = Data_mackerel_use_Ireland_select$julian_std,
                               y = Data_mackerel_use_Ireland_select$log_rate,
                               Likconfig = 0      # 0 = dnorm, 1 = dgamma
              )

              parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                                     log_sigma = rep(log(0.2),N_threshold)
              )

              Map = list()

              op <- getwd()
              setwd(paste0(getwd(),"/src"))
              obj2break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model_2breaks", map=Map)
              setwd(op)


              library(TMBhelper)
              opt_2breaks <- fit_tmb( obj=obj2break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
              opt_2breaks

              sd_report_2break <- sdreport(obj2break)
              check_estimability(obj2break)
              sigma <- as.vector(exp(summary(sd_report_2break, "fixed")[grep("log_sigma", rownames(summary(sd_report_2break, "fixed"))),1]))
              mu_pred <- matrix(summary(sd_report_2break, "report")[,1], ncol=3, byrow=FALSE)

              # Calculating the actual prediction
              LL <- obj2break$report()$LL


              # The weighting factor is the likelihood (but with the sum to 1 constraint)
              weight_2break <- LL
              weight_2break <- replace(weight_2break, which(weight_2break==0), NA)
              for (i in 1:data_tmb$N){
                weight_2break[i,,] <- exp(weight_2break[i,,])/sum(exp(weight_2break[i,,]), na.rm=T)
              }

              # Marginal 1st breakpoint value
              par(mfrow=c(2,1), mar=c(4,3,1,1), oma=c(1,1,1,1))
              plot(data_tmb$thresh, apply(apply(weight_2break,c(1,2),mean, na.rm=T), 2, mean), main ="First breakpoint", type="b")
              plot(data_tmb$thresh, apply(apply(weight_2break,c(1,3),mean, na.rm=T), 2, mean), main ="Second breakpoint", type="b")

              Prediction_2breaks <- rep(0, data_tmb$N)
              for (n in 1:data_tmb$N){
                for (thr in 1:(data_tmb$Nthres-1)) {
                  for (thr2 in (thr+1):data_tmb$Nthres) {
                    if (data_tmb$thres_cov[n] < (data_tmb$thresh[thr])){
                      Prediction_2breaks[n] = Prediction_2breaks[n]+weight_2break[n,thr,thr2]*mu_pred[n,1]
                    }
                    if ((data_tmb$thres_cov[n] >= (data_tmb$thresh[thr])) & (data_tmb$thres_cov[n] < (data_tmb$thresh[thr2]))){
                      Prediction_2breaks[n] = Prediction_2breaks[n]+weight_2break[n,thr,thr2]*mu_pred[n,2]
                    }
                    if (data_tmb$thres_cov[n] >= (data_tmb$thresh[thr2])){
                      Prediction_2breaks[n] = Prediction_2breaks[n]+weight_2break[n,thr,thr2]*mu_pred[n,3]
                    }
                  }
                }
              }

              plot(Prediction_2breaks, data_tmb$y); abline(0,1)
              qqnorm(y=(Prediction_2breaks-data_tmb$y)/sd(Prediction_2breaks-data_tmb$y))
              abline(0,1, lty=2)

              ## Now producing the residual vs predictor
              opt_2breaks_resid <- (Prediction_2breaks-data_tmb$y)/sd(Prediction_2breaks-data_tmb$y)
              par(mfrow=c(3,2))
              plot(Data_mackerel_use_Ireland_select$Length, opt_2breaks_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Tag_area), opt_2breaks_resid); abline(h=0, lty=2)
              plot(as.factor(Data_mackerel_use_Ireland_select$Release_timing), opt_2breaks_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$Release_year, opt_2breaks_resid); abline(h=0, lty=2)
              plot(Data_mackerel_use_Ireland_select$julian_recapture_std, opt_2breaks_resid); abline(h=0, lty=2)


              ## main effect table & plot
              Estimates_2break <- summary(sd_report_2break, "fixed")
              rownames(Estimates_2break) <- c(apply(cbind(colnames(XX),"before"),1, function(x) paste(x, collapse="_")),
                                              apply(cbind(colnames(XX),"mid"),1, function(x) paste(x, collapse="_")),
                                              apply(cbind(colnames(XX),"after"),1, function(x) paste(x, collapse="_")),
                                              "log_sigma_before", "log_sigma_mid", "log_sigma_after")
              Estimates_2break <- as.data.frame(Estimates_2break)
              Estimates_2break$gradient <- opt_2breaks$diagnostics$final_gradient

            }


          ## Simulating observations
				# 	Nit = 10000
				# 	Prediction <- array(0, dim=c(data_tmb$N, data_tmb$Nthres, Nit))
				# 	for (n in 1:data_tmb$N){
				# 	   for (thr in 1:data_tmb$Nthres){
				# 		  if (data_tmb$thres_cov[n] < (data_tmb$thresh[thr])){
				# 		  	Prediction[n,thr,] = rnorm(Nit, mu_pred[n,1], sigma[1])
				# 		  }
				# 		  if (data_tmb$thres_cov[n] >= (data_tmb$thresh[thr])){
				# 			  Prediction[n,thr,] = rnorm(Nit, mu_pred[n,2], sigma[2])
				# 		  }
				# 	   }
				# 	}
				#
				# 	Pred_val <- apply(Prediction, c(1,3), function(x) sum(weight*x))
				# 	Predict_val <- apply(Pred_val, 1, mean)
				# 	plot(Predict_val, data_tmb$y); abline(0,1)
				# 	qqnorm(y=(Predict_val-data_tmb$y)/sd(Predict_val-data_tmb$y))
				# 	abline(0,1, lty=2)
				#
				# # Conclusion:
				# Something weird is happening.... not the same results as the Stan model...
				# But it is technically the same model... just a different Latitudeform


      # grouping by year, area, and release timing
      dat_release_grouped  <- Data_mackerel_all %>% filter(Tag_area %in% c("South_Ireland", "North_Ireland"),
                                                  Release_year %in% 2014:2020) %>%
        mutate(Tag_area = droplevels(Tag_area), Release_year = droplevels(Release_year)) %>%
        group_by(Release_year, Tag_area, Release_timing, .drop = FALSE) %>% summarize(n_release=n())
      dat_recap_grouped  <- Data_mackerel_use_Ireland_select %>% mutate(Tag_area = droplevels(Tag_area)) %>%
        group_by(Release_year, Tag_area, Release_timing, .drop = FALSE) %>%
        summarize(n_recap=n())

      dat_recap_grouped$n_release <- dat_release_grouped$n_release
      dat_recap_grouped$recap_rate <- dat_recap_grouped$n_recap/dat_release_grouped$n_release

      # grouping by year, area
      dat_release_grouped1  <- Data_mackerel_all %>% filter(Tag_area %in% c("South_Ireland", "North_Ireland"),
                                                  Release_year %in% 2014:2020) %>%
        mutate(Tag_area = droplevels(Tag_area), Release_year = droplevels(Release_year)) %>%
        group_by(Release_year, Tag_area, .drop = FALSE) %>% summarize(n_release=n())
      dat_recap_grouped1  <- Data_mackerel_use_Ireland_select %>% mutate(Tag_area = droplevels(Tag_area)) %>%
        group_by(Release_year, Tag_area, .drop = FALSE) %>%
        summarize(n_recap=n())

      dat_recap_grouped1$n_release <- dat_release_grouped1$n_release
      dat_recap_grouped1$recap_rate <- dat_recap_grouped1$n_recap/dat_release_grouped1$n_release

      # grouping by year
      dat_release <- Data_mackerel_all %>% filter(Tag_area %in% c("South_Ireland", "North_Ireland"),
                                                  Release_year %in% 2014:2020) %>%
        mutate(Tag_area = droplevels(Tag_area), Release_year = droplevels(Release_year)) %>%
        group_by(Release_year, .drop = FALSE) %>% summarize(n_release=n())
      dat_recap <- Data_mackerel_use_Ireland_select %>% mutate(Tag_area = droplevels(Tag_area)) %>%
        group_by(Release_year, .drop = FALSE) %>%
        summarize(n_recap=n())

      dat_recap$n_release <- dat_release$n_release
      dat_recap$recap_rate <- dat_recap$n_recap/dat_release$n_release


      write.csv(dat_recap, file=paste0("../plots/to_aril/summary_recapture_coarse.csv"))
      write.csv(dat_recap_grouped1, file=paste0("../plots/to_aril/summary_recapture_year_area.csv"))
      write.csv(dat_recap_grouped, file=paste0("../plots/to_aril/summary_recapture_group.csv"))


      # For Northern Ireland
      map_area <- get_stamenmap(c(left = -32, bottom = 51, right = 21, top = 72), maptype = "toner-lite", zoom=5)
      g1 <- ggmap(map_area) +
        geom_jitter(data=subset(Data_mackerel_use_Ireland_select, Tag_area == "North_Ireland"), aes(x=cLon, y=cLat,col=duration), pch=16) +
        scale_color_viridis_c() + facet_wrap(. ~ Release_year) +
        scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6))
      g1

      g2 <- ggmap(map_area) +
        geom_jitter(data=subset(Data_mackerel_use_Ireland_select, Tag_area == "South_Ireland"), aes(x=cLon, y=cLat,col=duration), pch=16) +
        scale_color_viridis_c() + facet_wrap(. ~ Release_year) +
        scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6))
      g2

      Data_mackerel_use_Ireland_select %>% group_by(Release_timing) %>% summarize(quantile(la_scaled, seq(0.1,0.9,0.1)))


############# In the section below are some "nice" plots to investigate the mackerel distribution pattern
############# Using the tagging data

############ Some nice plots to look at mark recapture pattern


		#### Ivestigate the effect of release timing on the pattern of mackerel distribution
		plot_timerelease_simple <- function(Year, lag=0)
		{
				  dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
				                                              Tag_area %in% c("South_Ireland", "North_Ireland"))
				  dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
				  dat_release$Tag_area <- as.character(dat_release$Tag_area)
				  dat_release$Tag_area <- as.factor(dat_release$Tag_area)
				  dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
				  dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
				                                                Tag_area %in% c("South_Ireland", "North_Ireland"))
				  dat_recap$julian <-  as.numeric(julian(dat_recap$RecaptureDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
				  dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
				  dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)

				  dat_recap_new <- as_tibble(dat_recap) %>% dplyr::select(cLon, cLat, julian) %>% count(round(cLon,1), round(cLat,1), round(julian,0))
				  colnames(dat_recap_new) <- c("cLon", "cLat", "julian_recapture", "n")

				  dat_release_new <- as_tibble(dat_release) %>% dplyr::select(Longitude, Latitude, julian) %>% count(round(Longitude,1), round(Latitude,1), round(julian,0))
          colnames(dat_release_new) <- c("Longitude", "Latitude", "julian_release", "n")

				  ncount_release <- dat_release %>% count()
				  ncount_release$Longituden = c(-32)
				  ncount_release$Latitudet = c(53)
				  ncount_release$Label = paste0("n=", ncount_release$n)
				  ncount_recap <- dat_recap %>% count()
				  ncount_recap$Longituden = c(-16)
				  ncount_recap$Latitudet = c(70)
				  ncount_recap$Label = paste0("n=", ncount_recap$n, " (rate=", round(ncount_recap$n/ncount_release$n,3)*100, "%)")
				  Select_catch <- subset(Catch_data, subset=c(catchdate ==Year+lag))
				  Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
				  Select_catch <- subset(Select_catch, cLat<72)
				  hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)
				  hull <- dat_recap %>% dplyr::select(cLon, cLat) %>%  dplyr::slice(chull(cLon, cLat))
				  gravity <- hull %>%  mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
				    dplyr::select(mean_lon, mean_lat) %>% dplyr::distinct()
				  hull_new <- as.data.frame(hull)
				  coordinates(hull_new) <- c("cLon", "cLat")
				  crs(hull_new) <- "+proj=longlat +datum=WGS84"
				  new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
				  hull_new <- spTransform(hull_new, new_proj)
				  hull_new <- as.data.frame(hull_new)
				  hull_new1 <- hull_new

				  area <- Polygon(hull_new[,c('cLon','cLat')])@area
				  g1 <- ggplot(Norway) + geom_sf() + geom_point(data=dat_release_new, aes(x=Longitude, y=Latitude,col=julian_release, size=n)) +
				    #geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat)) +
				    scale_color_viridis_c() +
				    scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6)) +
				    new_scale_color() +
				    geom_point(data=dat_recap_new, aes(x=cLon, y=cLat, size=n)) +
				    scale_color_viridis_c() +
				    #scale_size_continuous(range = c(0.3, 0.6)) +
				    geom_text(data=data.frame(x=-8,y=72,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
				    geom_text(data=data.frame(x=-26,y=55,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
				    ggtitle(Year) +
				    geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA) #+
				  g1 <- g1 + geom_text(data=ncount_release, aes(x=lon, y=lat, label=label), hjust=0) +
				    geom_text(data=ncount_recap, aes(x=lon, y=lat, label=label), hjust=0) + labs(x="Longitude", y="Latitude") +
				    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
				  g1

				  RES <- list()
				  RES$plot <- g1
				  RES$area <- area
				  return(RES)
			}

		Make_plot_timing_simple <- function(Years, lag){
					  SAVE <- list()
					  # areas <- c()
					  for (yr in seq_along(Years)){
					    SAVE[[yr]] <- plot_timerelease_simple(Years[[yr]], lag=lag)
					    # areas <- rbind(areas, SAVE[[yr]][[2]])
					  }
					  # areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
					  # colnames(areas.std) <- Years
					  #
					  # areas_df <- reshape2::melt(areas.std, value.name="Size_area")
					  # colnames(areas_df) <- c("Year", "Size_area")
					  # areas_df$Year <- as.factor(areas_df$Year)
					  # areas_df$Time_release <- as.factor(areas_df$Time_release)
					  #
					  # lm1 <- (lm(Size_area ~ Time_release + Year, areas_df))

					  ## Plotting
					  ggg_timerelease_size <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
					  ggsave(ggg_timerelease_size, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)
					  # if(length(Years)<=4){
					  #   ggg_timerelease_size <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
					  #   ggsave(ggg_timerelease_size, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)
					  # }
					  # if(length(Years)>4 & length(Years)<=8){
					  #   ggg_timerelease_size1 <- arrangeGrob(grobs=map(SAVE[1:4], pluck, "plot"), nrow=2, ncol=2)
					  #   ggg_timerelease_size2 <- arrangeGrob(grobs=map(SAVE[5:length(Years)], pluck, "plot"), nrow=2, ncol=2)
					  #   ggsave(ggg_timerelease_size1, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, "_part1.pdf"), width=32, height=26, units="cm", device = "pdf", dpi = 400)
					  #   ggsave(ggg_timerelease_size2, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, "_part2.pdf"), width=32, height=26, units="cm", device = "pdf", dpi = 400)
					  # }
					  # if(length(Years)>8){
					  #   ggg_timerelease_size1 <- arrangeGrob(grobs=map(SAVE[1:4], pluck, "plot"), nrow=2, ncol=2)
					  #   ggg_timerelease_size2 <- arrangeGrob(grobs=map(SAVE[5:8], pluck, "plot"), nrow=2, ncol=2)
					  #   ggg_timerelease_size3 <- arrangeGrob(grobs=map(SAVE[9:length(Years)], pluck, "plot"), nrow=2, ncol=2)
					  #   ggsave(ggg_timerelease_size1, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, "_part1.pdf"), width=32, height=26, units="cm", device = "pdf", dpi = 400)
					  #   ggsave(ggg_timerelease_size2, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, "_part2.pdf"), width=32, height=26, units="cm", device = "pdf", dpi = 400)
					  #   ggsave(ggg_timerelease_size3, filename=paste0(getwd(), "/plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag", lag, "_part3.pdf"), width=32, height=26, units="cm", device = "pdf", dpi = 400)
					  # }

					  #RES <- list()
					  #RES$areas <- areas.std
					  #RES$areas.lm <- lm1
					  #return(RES)
					}

		# no lag : same year to 2 years lag in recapture
		(Area_timing_lag0 <- Make_plot_timing_simple(Years=2014:2020, lag=0))
		(Area_timing_lag1 <- Make_plot_timing_simple(Years=2011:2018, lag=1))
		(Area_timing_lag2 <- Make_plot_timing_simple(Years=2011:2017, lag=2))

		summary(Area_timing_lag0[[2]])
		summary(Area_timing_lag1[[2]])
		summary(Area_timing_lag2[[2]])





#### Ivestigate the effect of release region on the pattern of mackerel distribution
    plot_region <- function(Year, lag=0)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_recap$julian <-  as.numeric(julian(dat_recap$RecaptureDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))

      dat_recap_new <- as_tibble(dat_recap) %>% dplyr::select(cLon, cLat, Tag_area, julian) %>% count(round(cLon,1), round(cLat,1), Tag_area, round(julian,0))
      colnames(dat_recap_new) <- c("cLon", "cLat", "Tag_area", "julian_recapture", "n")

      dat_release_new <- as_tibble(dat_release) %>% dplyr::select(Longitude, Latitude, Tag_area, julian) %>% count(round(Longitude,1), round(Latitude,1), Tag_area, round(julian,0))
      colnames(dat_release_new) <- c("Longitude", "Latitude", "Tag_area", "julian_release", "n")

      ncount_release <- dat_release %>% group_by(Tag_area) %>% count()
      ncount_release$Longituden = c(-30, -30)
      ncount_release$Latitudet = c(58, 56.5)
      ncount_release$Label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(Tag_area) %>% count()
      ncount_recap$Longituden = c(-10, -10)
      ncount_recap$Latitudet = c(70, 69)
      ncount_recap$Label = paste0("n=", ncount_recap$n)

      dat_recap$Category <- paste(dat_recap$Tag_area, sep="_")

      Select_catch <- subset(Catch_data, subset=c(catchdate ==Year+lag))
      Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
      Select_catch <- subset(Select_catch, cLat<72)
      hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)

      hull <- dat_recap %>% dplyr::select(Tag_area, cLon, cLat) %>% group_by(Tag_area) %>%  dplyr::slice(chull(cLon, cLat))
      hull$Category <- paste(hull$Tag_area, sep="_")
      hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()

      cats <- c("From N. Ireland", "From S. Ireland")
      gravity <- hull %>% group_by(Category) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
      hull_new <- as.data.frame(hull)
      coordinates(hull_new) <- c("cLon", "cLat")
      crs(hull_new) <- "+proj=longlat +datum=WGS84"
      new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
      hull_new <- spTransform(hull_new, new_proj)
      hull_new <- as.data.frame(hull_new)
      Vertex_nb <- table(hull_new$Category)
      to_add <- which(Vertex_nb==3)
      to_remove <- which(Vertex_nb<3)
      if (length(to_add)>0){
        for (i in seq_along(to_add)){
          test <- subset(hull_new, Category==unique(hull_new$Category)[to_add[i]])
          test[4,] <- test[1,]
          hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_add[i]])
          hull_new <- rbind(hull_new, test)
        }
      }
      if (length(to_remove)>0){
        for (i in seq_along(to_remove)){
          hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_remove[i]])
        }
      }
      area <- sapply(unique(hull_new$Category), function(x) Polygon(subset(hull_new, Category==x)[,c('cLon','cLat')])@area)
      area_category <- rep(NA, length(cats))
      names(area_category) <- cats
      area_category[match(names(area), names(area_category))] <- area
      names(area_category) <- apply(cbind(cats, paste0("(n=", ncount_recap$n, ")")), 1, function(x) paste(x, collapse=" "))

#
#       prop_size <- dat_recap %>% group_by(Category,length_bin) %>% summarize(n=n()) %>% mutate(prop = prop.table(n)) %>% filter(n>1)
#       prop_size$r <- 1
#       Ngroups <- length(unique(prop_size$Category))
#       prop_size$Longituden <- -33
#       prop_size$Latitudet <- as.factor(prop_size$Category)
#       prop_size$Latitudet <- as.numeric(as.character(factor(prop_size$Latitudet, labels=seq(58,56+2.5*length(levels(prop_size$Latitudet)),by=2.5))))
#       prop_size_wide <- spread(prop_size[,-3], length_bin, prop)
#       prop_size_wide[is.na(prop_size_wide)] <- 0
#
#       prop_size_text <- data.frame(lon=-31, lat=seq(58,56+2.5*length(unique(prop_size$Latitudet)),by=2.5), label=unique(prop_size$Category))

      g1 <- ggmap(map_area) + geom_point(data=dat_release_new, aes(x=Longitude, y=Latitude,col=julian_release, size=n)) +
        scale_color_viridis_c() +
        scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6)) +
        new_scale_color() +
        geom_point(data=dat_recap_new, aes(x=cLon, y=cLat, size=n)) +
        scale_color_viridis_c() +
        geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set1")[c(1,2)]) +
        geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        geom_text(data=data.frame(x=-26,y=60,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) #+
        #geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA)
        g1 <- g1 + geom_text(data=ncount_release, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(3, "Set1")[c(1,2)]) +
          geom_text(data=ncount_recap, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(3, "Set1")[c(1,2)])

      RES <- list()
      RES$plot <- g1
      RES$area <- area_category
      return(RES)
    }

    Make_plot_region <- function(Years, lag){
      SAVE <- list()
      areas <- c()
      for (yr in seq_along(Years)){
        SAVE[[yr]] <- plot_region(Years[[yr]], lag=lag)
        areas <- rbind(areas, SAVE[[yr]][[2]])
      }
      areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
      colnames(areas.std) <- Years

      areas_df <- reshape2::melt(areas.std, value.name="Size_area")
      colnames(areas_df) <- c("Region", "Year", "Size_area")
      areas_df$Year <- as.factor(areas_df$Year)
      areas_df$Region <- as.factor(areas_df$Region)

      lm1 <- (lm(Size_area ~ Region + Year, areas_df))

      ## Plotting
      ggg_timerelease <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
      ggsave(ggg_timerelease, filename=paste0("../plots/mark_recap_region_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

      RES <- list()
      RES$areas <- areas.std
      RES$areas.lm <- lm1
      return(RES)
    }

    # no lag : same year to 2 years lag in recapture
    (Area_region_lag0 <- Make_plot_region(Years=2014:2020, lag=0))
    (Area_region_lag1 <- Make_plot_region(Years=2014:2018, lag=1))
    (Area_region_lag2 <- Make_plot_region(Years=2014:2017, lag=2))


#### Ivestigate the effect of release timing on the pattern of mackerel distribution
    plot_timerelease <- function(Year, lag=0)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))

      ncount_release <- dat_release %>% group_by(Release_timing) %>% count()
      ncount_release$Longituden = c(-30, -30)
      ncount_release$Latitudet = c(58, 56.5)
      ncount_release$Label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(Release_timing) %>% count()
      ncount_recap$Longituden = c(-10,-10)
      ncount_recap$Latitudet = c(70, 69)
      ncount_recap$Label = paste0("n=", ncount_recap$n)

      dat_recap$Category <- paste(dat_recap$Release_timing, sep="_")

      Select_catch <- subset(Catch_data, subset=c(catchdate ==Year+lag))
      Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
      Select_catch <- subset(Select_catch, cLat<72)
      hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)

      hull <- dat_recap %>% dplyr::select(Release_timing, cLon, cLat) %>% group_by(Release_timing) %>%  dplyr::slice(chull(cLon, cLat))
      hull$Category <- paste(hull$Release_timing, sep="_")
      hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
      #hull$Category <- as.factor(hull$Category)
      #hull$Category <- factor(hull$Category, levels=c("From N. Ireland_Early", "From N. Ireland_Mid",
      #                                                "From N. Ireland_Late", "From S. Ireland_Early",
      #                                                "From S. Ireland_Mid", "From S. Ireland_Late"))
      cats <- c("First_half", "Second_half")
      gravity <- hull %>% group_by(Category) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
      hull_new <- as.data.frame(hull)
      coordinates(hull_new) <- c("cLon", "cLat")
      crs(hull_new) <- "+proj=longlat +datum=WGS84"
      new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
      hull_new <- spTransform(hull_new, new_proj)
      hull_new <- as.data.frame(hull_new)
      hull_new1 <- hull_new
      Vertex_nb <- table(hull_new$Category)
      to_add <- which(Vertex_nb==3)
      to_remove <- which(Vertex_nb<3)
      if (length(to_add)>0){
        for (i in seq_along(to_add)){
          test <- subset(hull_new1, Category==unique(hull_new1$Category)[to_add[i]])
          test[4,] <- test[1,]
          hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_add[i]])
          hull_new <- rbind(hull_new, test)
        }
      }
      if (length(to_remove)>0){
        for (i in seq_along(to_remove)){
          hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_remove[i]])
        }
      }
      area <- sapply(unique(hull_new$Category), function(x) Polygon(subset(hull_new, Category==x)[,c('cLon','cLat')])@area)
      area_category <- rep(NA, length(cats))
      names(area_category) <- cats
      area_category[match(names(area), names(area_category))] <- area
      names(area_category) <- apply(cbind(cats, paste0("(n=", ncount_recap$n, ")")), 1, function(x) paste(x, collapse=" "))

      # prop_size <- dat_recap %>% group_by(Category,length_bin) %>% summarize(n=n()) %>% mutate(prop = prop.table(n)) %>% filter(n>1)
      # prop_size$r <- 1
      # Ngroups <- length(unique(prop_size$Category))
      # prop_size$Longituden <- -33
      # prop_size$Latitudet <- as.factor(prop_size$Category)
      # prop_size$Latitudet <- as.numeric(as.character(factor(prop_size$Latitudet, labels=seq(58,56+2.5*length(levels(prop_size$Latitudet)),by=2.5))))
      # prop_size_wide <- spread(prop_size[,-3], length_bin, prop)
      # prop_size_wide[is.na(prop_size_wide)] <- 0
      #
      # prop_size_text <- data.frame(lon=-31, lat=seq(58,56+2.5*length(unique(prop_size$Latitudet)),by=2.5), label=unique(prop_size$Category))

      g1 <- ggmap(map_area) + geom_jitter(data=dat_release, aes(x=Longitude, y=Latitude,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_point(data=dat_recap, aes(x=cLon, y=cLat), size=0.3, col="black") +
        geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set1")[c(1,2)]) +
        scale_color_viridis_c() +
        geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        geom_text(data=data.frame(x=-26,y=60,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) +
        geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA) #+
      # geom_scatterpie(aes(x=lon, y=lat, group = Category, r=r),
      #                data = prop_size_wide, cols = colnames(prop_size_wide[,-c(1:4)]))+
      # geom_text(data=prop_size_text, aes(x=lon, y=lat, label=label), size=2, hjust = 0)
      g1 <- g1 + geom_text(data=ncount_release, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(3, "Set1")[c(1,2)]) +
        geom_text(data=ncount_recap, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(3, "Set1")[c(1,2)])

      RES <- list()
      RES$plot <- g1
      RES$area <- area_category
      return(RES)
    }

    Make_plot_timing <- function(Years, lag){
      SAVE <- list()
      areas <- c()
      for (yr in seq_along(Years)){
        SAVE[[yr]] <- plot_timerelease(Years[[yr]], lag=lag)
        areas <- rbind(areas, SAVE[[yr]][[2]])
      }
      areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
      colnames(areas.std) <- Years

      areas_df <- reshape2::melt(areas.std, value.name="Size_area")
      colnames(areas_df) <- c("Time_release", "Year", "Size_area")
      areas_df$Year <- as.factor(areas_df$Year)
      areas_df$Time_release <- as.factor(areas_df$Time_release)

      lm1 <- (lm(Size_area ~ Time_release + Year, areas_df))

      ## Plotting
      ggg_timerelease_size <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
      ggsave(ggg_timerelease_size, filename=paste0("../plots/mark_recap_timing_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

      RES <- list()
      RES$areas <- areas.std
      RES$areas.lm <- lm1
      return(RES)
    }

    # no lag : same year to 2 years lag in recapture
      (Area_timing_lag0 <- Make_plot_timing(Years=2014:2020, lag=0))
      (Area_timing_lag1 <- Make_plot_timing(Years=2014:2018, lag=1))
      (Area_timing_lag2 <- Make_plot_timing(Years=2014:2017, lag=2))

      summary(Area_timing_lag0[[2]])
      summary(Area_timing_lag1[[2]])
      summary(Area_timing_lag2[[2]])


#### Ivestigate the effect of fish size on the pattern of mackerel distribution
    plot_fishsize <- function(Year, lag=0)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))

      ncount_release <- dat_release %>% group_by(length_bin) %>% count()
      ncount_release$Longituden = c(-32, -32, -32)
      ncount_release$Latitudet = c(58, 56.5, 55)
      ncount_release$Label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(length_bin) %>% count()
      ncount_recap$Longituden = c(-10, -10,-10)
      ncount_recap$Latitudet = c(70, 69, 68)
      ncount_recap$Label = paste0("n=", ncount_recap$n)

      dat_recap$Category <- paste(dat_recap$length_bin, sep="_")

      Select_catch <- subset(Catch_data, subset=c(catchdate ==Year+lag))
      Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
      Select_catch <- subset(Select_catch, cLat<72)
      hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)

      hull <- dat_recap %>% dplyr::select(length_bin, cLon, cLat) %>% group_by(length_bin) %>%  dplyr::slice(chull(cLon, cLat))
      hull$Category <- paste(hull$length_bin, sep="_")
      hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
      #hull$Category <- as.factor(hull$Category)
      #hull$Category <- factor(hull$Category, levels=c("From N. Ireland_Early", "From N. Ireland_Mid",
      #                                                "From N. Ireland_Late", "From S. Ireland_Early",
      #                                                "From S. Ireland_Mid", "From S. Ireland_Late"))
      cats <- c("(0,34cm]", "(34cm,37cm]","(37cm,45cm]")
      gravity <- hull %>% group_by(Category) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
      hull_new <- as.data.frame(hull)
      coordinates(hull_new) <- c("cLon", "cLat")
      crs(hull_new) <- "+proj=longlat +datum=WGS84"
      new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
      hull_new <- spTransform(hull_new, new_proj)
      hull_new <- as.data.frame(hull_new)
      hull_new1 <- hull_new
      Vertex_nb <- table(hull_new$Category)
      to_add <- which(Vertex_nb==3)
      to_remove <- which(Vertex_nb<3)
      if (length(to_add)>0){
        for (i in seq_along(to_add)){
          test <- subset(hull_new1, Category==unique(hull_new1$Category)[to_add[i]])
          test[4,] <- test[1,]
          hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_add[i]])
          hull_new <- rbind(hull_new, test)
        }
      }
      if (length(to_remove)>0){
        for (i in seq_along(to_remove)){
          hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_remove[i]])
        }
      }
      area <- sapply(unique(hull_new$Category), function(x) Polygon(subset(hull_new, Category==x)[,c('cLon','cLat')])@area)
      area_category <- rep(NA, length(cats))
      names(area_category) <- cats
      area_category[match(names(area), names(area_category))] <- area
      names(area_category) <- apply(cbind(cats, paste0("(n=", ncount_recap$n, ")")), 1, function(x) paste(x, collapse=" "))

      # prop_size <- dat_recap %>% group_by(Category,length_bin) %>% summarize(n=n()) %>% mutate(prop = prop.table(n)) %>% filter(n>1)
      # prop_size$r <- 1
      # Ngroups <- length(unique(prop_size$Category))
      # prop_size$Longituden <- -33
      # prop_size$Latitudet <- as.factor(prop_size$Category)
      # prop_size$Latitudet <- as.numeric(as.character(factor(prop_size$Latitudet, labels=seq(58,56+2.5*length(levels(prop_size$Latitudet)),by=2.5))))
      # prop_size_wide <- spread(prop_size[,-3], length_bin, prop)
      # prop_size_wide[is.na(prop_size_wide)] <- 0
      #
      # prop_size_text <- data.frame(lon=-31, lat=seq(58,56+2.5*length(unique(prop_size$Latitudet)),by=2.5), label=unique(prop_size$Category))

      g1 <- ggmap(map_area) + geom_jitter(data=dat_release, aes(x=Longitude, y=Latitude,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_point(data=dat_recap, aes(x=cLon, y=cLat), size=0.3, col="black") +
        geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set1")[1:3]) +
        scale_color_viridis_c() +
        geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        geom_text(data=data.frame(x=-25,y=60,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) +
        geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA) #+
      # geom_scatterpie(aes(x=lon, y=lat, group = Category, r=r),
      #                data = prop_size_wide, cols = colnames(prop_size_wide[,-c(1:4)]))+
      # geom_text(data=prop_size_text, aes(x=lon, y=lat, label=label), size=2, hjust = 0)
        g1 <- g1 + geom_text(data=ncount_release, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(3, "Set1")[1:3]) +
        geom_text(data=ncount_recap, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(3, "Set1")[1:3])

      RES <- list()
      RES$plot <- g1
      RES$area <- area_category
      return(RES)
    }

    Make_plot_fishsize <- function(Years, lag){
      SAVE <- list()
      areas <- c()
      for (yr in seq_along(Years)){
        SAVE[[yr]] <- plot_fishsize(Years[[yr]], lag=lag)
        areas <- rbind(areas, SAVE[[yr]][[2]])
      }
      areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
      colnames(areas.std) <- Years

      areas_df <- reshape2::melt(areas.std, value.name="Size_area")
      colnames(areas_df) <- c("Fish_size", "Year", "Size_area")
      areas_df$Year <- as.factor(areas_df$Year)
      areas_df$Fish_size <- as.factor(areas_df$Fish_size)

      lm1 <- (lm(Size_area ~ Fish_size + Year, areas_df))

      ## Plotting
      ggg_timerelease_size <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
      ggsave(ggg_timerelease_size, filename=paste0("../plots/mark_recap_size_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

      RES <- list()
      RES$areas <- areas.std
      RES$areas.lm <- lm1
      return(RES)
    }

    # no lag : same year to 2 years lag in recapture
      (Area_fishsize_lag0 <- Make_plot_fishsize(Years=2014:2020, lag=0))
      (Area_fishsize_lag1 <- Make_plot_fishsize(Years=2014:2018, lag=1))
      (Area_fishsize_lag2 <- Make_plot_fishsize(Years=2014:2017, lag=2))

      summary(Area_fishsize_lag0[[2]])
      summary(Area_fishsize_lag1[[2]])
      summary(Area_fishsize_lag2[[2]])


#### Plot to investigate the location of recapture from Ireland, by year, time of release, region and fish size
     plot_region_timerelease <- function(Year, lag=0)
     {
        dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                      Tag_area %in% c("South_Ireland", "North_Ireland"))
        dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
        dat_release$Tag_area <- as.character(dat_release$Tag_area)
        dat_release$Tag_area <- as.factor(dat_release$Tag_area)
        dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
        dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                                      Tag_area %in% c("South_Ireland", "North_Ireland"))
        dat_recap$julian <-  as.numeric(julian(dat_recap$RecaptureDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
        dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
        dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
        dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
        # dat_recap <- dat_recap %>% mutate(Release_timing = ifelse(Release_monthday <= "5_14", "Early",
        # ifelse((Release_monthday <= "5_31" & Release_monthday >= "5_15"), "Mid", "Late")))

        dat_release$Category <- paste(dat_release$Tag_area, dat_release$Release_timing, sep="_")
        dat_recap$Category <- paste(dat_recap$Tag_area, dat_recap$Release_timing, sep="_")

        dat_recap_new <- as_tibble(dat_recap) %>% dplyr::select(cLon, cLat, Tag_area, Release_timing, julian) %>% count(round(cLon,1), round(cLat,1), Tag_area, Release_timing, round(julian,0))
        colnames(dat_recap_new) <- c("cLon", "cLat", "Tag_area", "Release_timing", "julian_recapture", "n")

        dat_release_new <- as_tibble(dat_release) %>% dplyr::select(Longitude, Latitude, Tag_area, Release_timing, julian) %>% count(round(Longitude,1), round(Latitude,1), Tag_area, Release_timing, round(julian,0))
        colnames(dat_release_new) <- c("Longitude", "Latitude", "Tag_area", "Release_timing", "julian_release", "n")

        ncount_recap <- dat_recap %>% group_by(Category) %>% count()
        morethan3 <- which(ncount_recap$n >= 3)
        ncount_recap <- ncount_recap[morethan3,]
        ncount_recap$Longituden = rep(-10,nrow(ncount_recap))
        ncount_recap$Latitudet = seq(70,(70-(nrow(ncount_recap)-1)), by=-1)
        ncount_recap$Label = paste0("n=", ncount_recap$n)

        ncount_release <- dat_release %>% group_by(Category) %>% count()
        morethan3 <- match(ncount_recap$Category, ncount_release$Category)
        ncount_release <- ncount_release[morethan3,]
        ncount_release$Longituden = rep(-32,nrow(ncount_release))
        ncount_release$Latitudet = seq(58,(58-1.5*(nrow(ncount_release)-1)), by=-1.5)
        ncount_release$Label = paste0("n=", ncount_release$n)

        Select_catch <- subset(Catch_data, subset=c(catchdate ==Year+lag))
        Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
        Select_catch <- subset(Select_catch, cLat<72)
        hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)

        hull <- dat_recap %>% dplyr::select(Tag_area, Release_timing, cLon, cLat) %>% group_by(Tag_area, Release_timing) %>%  dplyr::slice(chull(cLon, cLat))
        hull$Category <- paste(hull$Tag_area, hull$Release_timing, sep="_")
        hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
        #hull$Category <- as.factor(hull$Category)
        #hull$Category <- factor(hull$Category, levels=c("From N. Ireland_Early", "From N. Ireland_Mid",
        #                                                "From N. Ireland_Late", "From S. Ireland_Early",
        #                                                "From S. Ireland_Mid", "From S. Ireland_Late"))
        cats <- c("From N. Ireland_First_half", "From N. Ireland_Second_half",
                  "From S. Ireland_First_half", "From S. Ireland_Second_half")
        gravity <- hull %>% group_by(Category) %>%
          mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
          dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
        hull_new <- as.data.frame(hull)
        coordinates(hull_new) <- c("cLon", "cLat")
        crs(hull_new) <- "+proj=longlat +datum=WGS84"
        new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
        hull_new <- spTransform(hull_new, new_proj)
        hull_new <- as.data.frame(hull_new)
        hull_new1 <- hull_new
        Vertex_nb <- table(hull_new$Category)
        to_add <- which(Vertex_nb==3)
        to_remove <- which(Vertex_nb<3)
        if (length(to_add)>0){
          for (i in seq_along(to_add)){
            test <- subset(hull_new1, Category==unique(hull_new1$Category)[to_add[i]])
            test[4,] <- test[1,]
            hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_add[i]])
            hull_new <- rbind(hull_new, test)
          }
        }
        if (length(to_remove)>0){
          for (i in seq_along(to_remove)){
            hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_remove[i]])
          }
        }
        area <- sapply(unique(hull_new$Category), function(x) Polygon(subset(hull_new, Category==x)[,c('cLon','cLat')])@area)
        area_category <- rep(NA, length(cats))
        names(area_category) <- cats
        area_category[match(names(area), names(area_category))] <- area

        # prop_size <- dat_recap %>% group_by(Category,length_bin) %>% summarize(n=n()) %>% mutate(prop = prop.table(n)) %>% filter(n>1)
        # prop_size$r <- 1
        # Ngroups <- length(unique(prop_size$Category))
        # prop_size$Longituden <- -33
        # prop_size$Latitudet <- as.factor(prop_size$Category)
        # prop_size$Latitudet <- as.numeric(as.character(factor(prop_size$Latitudet, labels=seq(58,56+2.5*length(levels(prop_size$Latitudet)),by=2.5))))
        # prop_size_wide <- spread(prop_size[,-3], length_bin, prop)
        # prop_size_wide[is.na(prop_size_wide)] <- 0
        #
        # prop_size_text <- data.frame(lon=-31, lat=seq(58,56+2.5*length(unique(prop_size$Latitudet)),by=2.5), label=unique(prop_size$Category))

        g1 <- ggmap(map_area) + geom_point(data=dat_release_new, aes(x=Longitude, y=Latitude,col=julian_release, size=n)) +
          scale_color_viridis_c() +
          scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6)) +
          new_scale_color() +
          geom_point(data=dat_recap_new, aes(x=cLon, y=cLat, size=n)) +
          scale_color_viridis_c() +
          geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
          scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Set1")[1:nrow(ncount_release)]) +
          geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
          geom_text(data=data.frame(x=-25,y=60,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
          ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) +
          scale_size(guide = 'none')
          #geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA) #+
        # geom_scatterpie(aes(x=lon, y=lat, group = Category, r=r),
        #                data = prop_size_wide, cols = colnames(prop_size_wide[,-c(1:4)]))+
        # geom_text(data=prop_size_text, aes(x=lon, y=lat, label=label), size=2, hjust = 0)
        g1 <- g1 + geom_text(data=ncount_release, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(4, "Set1")[1:nrow(ncount_release)]) +
          geom_text(data=ncount_recap, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(4, "Set1")[1:nrow(ncount_release)])

        RES <- list()
        RES$plot <- g1
        RES$area <- area_category
        return(RES)
      }

     Make_plot_region_timerelease <- function(Years, lag){
        SAVE <- list()
        areas <- c()
        for (yr in seq_along(Years)){
          SAVE[[yr]] <- plot_region_timerelease(Years[[yr]], lag=lag)
          areas <- rbind(areas, SAVE[[yr]][[2]])
        }
        areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
        colnames(areas.std) <- Years

        areas_df <- reshape2::melt(areas.std, value.name="Size_area")
        colnames(areas_df) <- c("Category", "Year", "Size_area")
        areas_df$Year <- as.factor(areas_df$Year)
        areas_df$Category <- as.factor(areas_df$Category)

        lm1 <- (lm(Size_area ~ Category + Year, areas_df))

        ## Plotting
        ggg_timerelease_size <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
        ggsave(ggg_timerelease_size, filename=paste0("../plots/mark_recap_region_timing_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

        RES <- list()
        RES$areas <- areas.std
        RES$areas.lm <- lm1
        return(RES)
      }
     # no lag : same year to 2 years lag in recapture
     (Area_region_timing_lag0 <- Make_plot_region_timerelease(Years=2014:2020, lag=0))
     (Area_region_timing_lag1 <- Make_plot_region_timerelease(Years=2014:2018, lag=1))
     (Area_region_timing_lag2 <- Make_plot_region_timerelease(Years=2014:2017, lag=2))
     summary(Area_region_timing_lag0[[2]])
     summary(Area_region_timing_lag1[[2]])
     summary(Area_region_timing_lag2[[2]])


#### Plot to investigate the location of recapture from Ireland, by year, time of release, region and fish size
     plot_region_time_size <- function(Year, lag=0)
     {
       dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                     Tag_area %in% c("South_Ireland", "North_Ireland"))
       dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
       dat_release$Tag_area <- as.character(dat_release$Tag_area)
       dat_release$Tag_area <- as.factor(dat_release$Tag_area)
       dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
       dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                                     Tag_area %in% c("South_Ireland", "North_Ireland"))
       dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
       dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
       dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
       # dat_recap <- dat_recap %>% mutate(Release_timing = ifelse(Release_monthday <= "5_14", "Early",
       # ifelse((Release_monthday <= "5_31" & Release_monthday >= "5_15"), "Mid", "Late")))

       dat_release$Category <- paste(dat_release$Tag_area, dat_release$Release_timing, dat_release$length_bin, sep="_")
       dat_recap$Category <- paste(dat_recap$Tag_area, dat_recap$Release_timing, dat_recap$length_bin, sep="_")

       ncount_recap <- dat_recap %>% group_by(Category) %>% count()
       morethan3 <- which(ncount_recap$n >= 3)
       ncount_recap <- ncount_recap[morethan3,]
       if (nrow(ncount_recap) <4) ncount_recap$Longituden = rep(-10,nrow(ncount_recap))
       if (nrow(ncount_recap) >=4 & nrow(ncount_recap) <8) ncount_recap$Longituden = c(rep(-13,4),rep(-5,nrow(ncount_recap)%%4))
       if (nrow(ncount_recap) >=8) ncount_recap$Longituden = c(rep(-14,4),rep(-6, 4), rep(2,nrow(ncount_recap)%%4))
       if (nrow(ncount_recap) <4) ncount_recap$Latitudet = seq(70,(70-(nrow(ncount_recap)-1)), by=-1)
       if (nrow(ncount_recap) >=4 & nrow(ncount_recap) <8) ncount_recap$Latitudet = c(70:67, seq(70,(70-(nrow(ncount_recap)%%4-1)), by=-1))
       if (nrow(ncount_recap) ==8) ncount_recap$Latitudet = c(70:67, 70:67)
       if (nrow(ncount_recap) >=9) ncount_recap$Latitudet = c(70:67, 70:67, seq(70,(70-(nrow(ncount_recap)%%4-1)), by=-1))
       ncount_recap$Label = paste0("n=", ncount_recap$n)

       ncount_release <- dat_release %>% group_by(Category) %>% count()
       morethan3 <- match(ncount_recap$Category, ncount_release$Category)
       ncount_release <- ncount_release[morethan3,]
       if (nrow(ncount_release) <5) ncount_release$Longituden = rep(-35,nrow(ncount_release))
       if (nrow(ncount_release) >=5 & nrow(ncount_release) <10) ncount_release$Longituden = c(rep(-35,5),rep(-24,nrow(ncount_release)%%5))
       if (nrow(ncount_release) >=10) ncount_release$Longituden = c(rep(-35,5),rep(-24,5), rep(-14,nrow(ncount_release)%%5))
       if (nrow(ncount_release) <5) ncount_release$Latitudet = seq(58,(58-1.5*(nrow(ncount_release)-1)), by=-1.5)
       if (nrow(ncount_release) >=5 & nrow(ncount_release) <10) ncount_release$Latitudet = c(seq(58,52,by=-1.5), seq(58,(58-1.5*(nrow(ncount_release)%%5-1)), by=-1.5))
       if (nrow(ncount_release) ==10) ncount_release$Latitudet = c(seq(58,52,by=-1.5), seq(58,52,by=-1.5))
       if (nrow(ncount_release) >=11) ncount_release$Latitudet = c(seq(58,52,by=-1.5), seq(58,52,by=-1.5), seq(58,(58-1.5*(nrow(ncount_release)%%5-1)), by=-1.5))
       ncount_release$Label = paste0("n=", ncount_release$n)

       Select_catch <- subset(Catch_data, subset=c(catchdate ==Year+lag))
       Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
       Select_catch <- subset(Select_catch, cLat<72)
       hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)

       hull <- dat_recap %>% dplyr::select(Tag_area, Release_timing, length_bin, cLon, cLat) %>% group_by(Tag_area, Release_timing, length_bin) %>%  dplyr::slice(chull(cLon, cLat))
       hull$Category <- paste(hull$Tag_area, hull$Release_timing, hull$length_bin, sep="_")
       hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
       #hull$Category <- as.factor(hull$Category)
       #hull$Category <- factor(hull$Category, levels=c("From N. Ireland_Early", "From N. Ireland_Mid",
       #                                                "From N. Ireland_Late", "From S. Ireland_Early",
       #                                                "From S. Ireland_Mid", "From S. Ireland_Late"))
       cats <- c("From N. Ireland_First_half", "From N. Ireland_Second_half",
                 "From S. Ireland_First_half", "From S. Ireland_Second_half")
       size <- levels(Data_mackerel_use_Ireland_select$length_bin)
       cats <- expand.grid(cats, size)
       cats <- apply(cats, 1, function(x) paste(x, collapse="_"))
       gravity <- hull %>% group_by(Category) %>%
         mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
         dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
       hull_new <- as.data.frame(hull)
       coordinates(hull_new) <- c("cLon", "cLat")
       crs(hull_new) <- "+proj=longlat +datum=WGS84"
       new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
       hull_new <- spTransform(hull_new, new_proj)
       hull_new <- as.data.frame(hull_new)
       hull_new1 <- hull_new
       Vertex_nb <- table(hull_new$Category)
       to_add <- which(Vertex_nb==3)
       to_remove <- which(Vertex_nb<3)
       if (length(to_add)>0){
         for (i in seq_along(to_add)){
           test <- subset(hull_new1, Category==unique(hull_new1$Category)[to_add[i]])
           test[4,] <- test[1,]
           hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_add[i]])
           hull_new <- rbind(hull_new, test)
         }
       }
       if (length(to_remove)>0){
         for (i in seq_along(to_remove)){
           hull_new <- subset(hull_new, Category!=unique(hull_new$Category)[to_remove[i]])
         }
       }
       area <- sapply(unique(hull_new$Category), function(x) Polygon(subset(hull_new, Category==x)[,c('cLon','cLat')])@area)
       area_category <- rep(NA, length(cats))
       names(area_category) <- cats
       area_category[match(names(area), names(area_category))] <- area

       # prop_size <- dat_recap %>% group_by(Category,length_bin) %>% summarize(n=n()) %>% mutate(prop = prop.table(n)) %>% filter(n>1)
       # prop_size$r <- 1
       # Ngroups <- length(unique(prop_size$Category))
       # prop_size$Longituden <- -33
       # prop_size$Latitudet <- as.factor(prop_size$Category)
       # prop_size$Latitudet <- as.numeric(as.character(factor(prop_size$Latitudet, labels=seq(58,56+2.5*length(levels(prop_size$Latitudet)),by=2.5))))
       # prop_size_wide <- spread(prop_size[,-3], length_bin, prop)
       # prop_size_wide[is.na(prop_size_wide)] <- 0
       #
       # prop_size_text <- data.frame(lon=-31, lat=seq(58,56+2.5*length(unique(prop_size$Latitudet)),by=2.5), label=unique(prop_size$Category))

       g1 <- ggmap(map_area) + geom_jitter(data=dat_release, aes(x=Longitude, y=Latitude,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
         geom_point(data=dat_recap, aes(x=cLon, y=cLat), size=0.3, col="black") +
         geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
         scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3")[1:nrow(ncount_release)]) +
         scale_color_viridis_c() +
         geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
         geom_text(data=data.frame(x=-25,y=60,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
         ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) +
         geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA) #+
       # geom_scatterpie(aes(x=lon, y=lat, group = Category, r=r),
       #                data = prop_size_wide, cols = colnames(prop_size_wide[,-c(1:4)]))+
       # geom_text(data=prop_size_text, aes(x=lon, y=lat, label=label), size=2, hjust = 0)
       g1 <- g1 + geom_text(data=ncount_release, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(12, "Set3")[1:nrow(ncount_release)]) +
         geom_text(data=ncount_recap, aes(x=lon, y=lat, label=label), hjust=0, col = RColorBrewer::brewer.pal(12, "Set3")[1:nrow(ncount_release)])

       RES <- list()
       RES$plot <- g1
       RES$area <- area_category
       return(RES)
     }

     Make_plot_region_time_size <- function(Years, lag){
       SAVE <- list()
       areas <- c()
       for (yr in seq_along(Years)){
         SAVE[[yr]] <- plot_region_time_size(Years[[yr]], lag=lag)
         areas <- rbind(areas, SAVE[[yr]][[2]])
       }
       areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
       colnames(areas.std) <- Years

       areas_df <- reshape2::melt(areas.std, value.name="Size_area")
       colnames(areas_df) <- c("Category", "Year", "Size_area")
       areas_df$Year <- as.factor(areas_df$Year)
       areas_df$Category <- as.factor(areas_df$Category)

       lm1 <- (lm(Size_area ~ Category + Year, areas_df))

       ## Plotting
       ggg_timerelease_size <- arrangeGrob(grobs=map(SAVE, pluck, "plot"), nrow=3, ncol=2)
       ggsave(ggg_timerelease_size, filename=paste0("../plots/mark_recap_region_timing_size_", Years[1], "_", tail(Years, 1), "_lag", lag, ".pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

       RES <- list()
       RES$areas <- areas.std
       RES$areas.lm <- lm1
       return(RES)
     }

     # no lag : same year to 2 years lag in recapture
     (Area_region_time_size_lag0 <- Make_plot_region_time_size(Years=2014:2020, lag=0))
     (Area_region_time_size_lag1 <- Make_plot_region_time_size(Years=2014:2018, lag=1))
     (Area_region_time_size_lag2 <- Make_plot_region_time_size(Years=2014:2017, lag=2))
     summary(Area_region_time_size_lag0[[2]])
     summary(Area_region_time_size_lag1[[2]])
     summary(Area_region_time_size_lag2[[2]])




  ### Plot to investigate the location of recapture from Ireland, by year and
  ### Duration between release and catch
  ### To Determine if what matter is time of release or duration
    plot_Yr_duration <- function(Year)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_release$julian <-  as.numeric(julian(dat_release$ReleaseDate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year &
                                                    Tag_area %in% c("South_Ireland", "North_Ireland"))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From S. Ireland"))

      dat_recap$duration_bin <- cut(dat_recap$duration, breaks=c(0,100,150,250))

      ncount_release <- dat_release %>% group_by(Tag_area) %>% count()
      ncount_release$Longituden = c(-17, -17)
      ncount_release$Latitudet = c(57, 51)
      ncount_release$Label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(Tag_area) %>% count()
      ncount_recap$Longituden = c(-8, -8)
      ncount_recap$Latitudet = c(69, 70)
      ncount_recap$Label = paste0("n=", ncount_recap$n)


      hull <- dat_recap %>% dplyr::select(Tag_area, duration_bin, cLon, cLat) %>% group_by(Tag_area, duration_bin) %>%  dplyr::slice(chull(cLon, cLat))
      hull$Category <- paste(hull$Tag_area, hull$duration_bin, sep="_")
      hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
      gravity <- hull %>% group_by(Category) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
      #area <- hull %>% group_by(Category) %>% Polygon()

      g1 <- ggmap(map_area) + geom_jitter(data=subset(dat_release, Tag_area=="From S. Ireland"), aes(x=Longitude, y=Latitude,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_point(data=subset(dat_recap, Tag_area=="From S. Ireland"), aes(x=cLon, y=cLat), size=0.3, col="black") +
        geom_jitter(data=subset(dat_release, Tag_area=="From N. Ireland"), aes(x=Longitude, y=Latitude,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_text(data=subset(ncount_recap, Tag_area=="From S. Ireland"), aes(x=lon, y=lat, label=label), col="red") +
        geom_text(data=subset(ncount_release, Tag_area=="From S. Ireland"), aes(x=lon, y=lat, label=label), col="red")  +
        geom_text(data=subset(ncount_recap, Tag_area=="From N. Ireland"), aes(x=lon, y=lat, label=label), col="blue") +
        geom_text(data=subset(ncount_release, Tag_area=="From N. Ireland"), aes(x=lon, y=lat, label=label), col="blue") +
        geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
        scale_color_viridis_c() +
        geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        geom_text(data=data.frame(x=-25,y=54,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5))
      g1

      return(g1)
    }

    ## 2011-2013, there is no recapture originating from both South and North side of Ireland
    p2014_duration <- plot_Yr_duration(2014)
    p2015_duration <- plot_Yr_duration(2015)
    p2016_duration <- plot_Yr_duration(2016)
    p2017_duration <- plot_Yr_duration(2017)
    p2018_duration <- plot_Yr_duration(2018)
    p2019_duration <- plot_Yr_duration(2019)

    ggg_duration <- arrangeGrob(p2014_duration,p2015_duration,p2016_duration,
                                   p2017_duration,p2018_duration,p2019_duration,
                                   nrow=3, ncol=2)
    ggsave(ggg_duration, filename=paste0("../plots/mark_recap_duration_2014-2019.pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)





################# Some other analyses

  #### Limiting the analysis to releases around Iceland:

      Data_mackerel_use_Iceland <- subset(Data_mackerel_final, Tag_area %in% c("Iceland"))

      ## Using duration (not standardized) but focusing on within year recapture
      Data_mackerel_use_Iceland_select <- subset(Data_mackerel_use_Iceland, duration < 365)
      m0_Iceland <- gam(log_rate ~ length + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Iceland_select)
      m1_Iceland <- gam(log_rate ~ length + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_Iceland_select)
      m2_Iceland <- gam(log_rate ~ length + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Iceland_select)
      m3_Iceland <- gam(log_rate ~ length + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Iceland_select)
      m4_Iceland <- gam(log_rate ~ length + Release_year, family=gaussian, data=Data_mackerel_use_Iceland_select)
      m5_Iceland <- gam(log_rate ~ length, family=gaussian, data=Data_mackerel_use_Iceland_select)

      AICtab(m1_Iceland,m2_Iceland,m3_Iceland,m4_Iceland,m5_Iceland)
      plot(m4_Iceland, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_select_Iceland <- simulateResiduals(fittedModel = m4_Iceland, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_select_Iceland, rank=T, quantreg = TRUE)
      testResiduals(sim_select_Iceland)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Iceland_select$Length, sim_select_Iceland$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Iceland_select$Release_month, sim_select_Iceland$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Iceland_select$Release_year, sim_select_Iceland$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Iceland_select$Latitude, sim_select_Iceland$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Iceland_select$Longitude, sim_select_Iceland$scaledResiduals, main = "Longitude")

      ## Using duration std1
      mmm0_Iceland <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm1_Iceland <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude) + s(date), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm2_Iceland <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Iceland)
      mmm3_Iceland <- gam(log_rate_annual1 ~ s(length, by= Release_year) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm4_Iceland <- gam(log_rate_annual1 ~ s(Longitude, Latitude, by= Release_year) + s(length), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm5_Iceland <- gam(log_rate_annual1 ~ s(length) + Release_year, family=gaussian, data=Data_mackerel_use_Iceland)

      AICtab(mmm0_Iceland, mmm1_Iceland,mmm2_Iceland,mmm3_Iceland,mmm4_Iceland,mmm5_Iceland,mmm6_Iceland)
      plot(mmm3_Iceland, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_std1_Iceland <- simulateResiduals(fittedModel = mmm3_Iceland, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_std1_Iceland, rank=T, quantreg = TRUE)
      testResiduals(sim_std1_Iceland)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Iceland$Length, sim_std1_Iceland$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Iceland$Release_month, sim_std1_Iceland$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Iceland$Release_year, sim_std1_Iceland$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Iceland$Latitude, sim_std1_Iceland$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Iceland$Longitude, sim_std1_Iceland$scaledResiduals, main = "Longitude")



  #### Limiting the analysis to releases around Bergen:

      Data_mackerel_use_Bergen <- subset(Data_mackerel_final, Tag_area %in% c("Bergen"))

      ## Using duration std1
      mmm0_Bergen <- gam(log_rate_annual1 ~ 1, family=gaussian, data=Data_mackerel_use_Bergen)
      mmm1_Bergen <- gam(log_rate_annual1 ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Bergen)
      mmm2_Bergen <- gam(log_rate_annual1 ~ s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Bergen)
      mmm3_Bergen <- gam(log_rate_annual1 ~ s(length), family=gaussian, data=Data_mackerel_use_Bergen)

      AICtab(mmm0_Bergen, mmm1_Bergen,mmm2_Bergen,mmm3_Bergen)
      plot(mmm1_Bergen, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_std1_Bergen <- simulateResiduals(fittedModel = mmm1_Bergen, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_std1_Bergen, rank=T, quantreg = TRUE)
      testResiduals(sim_std1_Bergen)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Bergen$Length, sim_std1_Bergen$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Bergen$Release_month, sim_std1_Bergen$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Bergen$Release_year, sim_std1_Bergen$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Bergen$Latitude, sim_std1_Bergen$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Bergen$Longitude, sim_std1_Bergen$scaledResiduals, main = "Longitude")



#### Maybe some tag recapture modeling?
## What can a mark recapture model provide us?

# Petersen model = closed population model (so far)
# brownie model

## survival rate, total mortality (fishing+natural), movement probability between areas
## rfid = pit tag
## high total mortality (80-90%) --> which can be related to tag loss/shedding
## Tagging mortality rate = see above
## Tag reporting rate = see above
## No info on fish size at the time of capture? No
## Is it all the catch data? (not all the factories, and some do no scan for rfid, some factories might not be representative/issues --> Aril has some documents on it)

#      1. migration area too, capture month, weight how much scanned
#      2. icelanding start fishin in July-August
#      3. September-onwards norwegian (fish start to move south  - roughly are finished feeding migration north)

# --> Based on Subbey et al. report, this might not be very informative



#
#   #### Limiting the analysis to releases around Ireland:
#
#   Data_mackerel_use_Ireland <- subset(Data_mackerel_final, Tag_area %in% c("South_Ireland", "North_Ireland"))
#
#   ## Using duration (not standardized) but focusing on within year recapture
#   Data_mackerel_use_Ireland_select <- subset(Data_mackerel_use_Ireland, duration < 200)
#   Data_mackerel_use_Ireland_select$Longitudecation <- with(Data_mackerel_use_Ireland_select, glmmTMB::numFactor(Longitude,Latitude))
#   levels(Data_mackerel_use_Ireland_select$Longitudecation)
#   Data_mackerel_use_Ireland_select$group <- factor(rep(1, nrow(Data_mackerel_use_Ireland_select)))
#   # Movement East_South
#     library(glmmTMB)
#     #mm0 <- glmmTMB(EW_move/1000 ~ offset(duration) + length + Tag_area, family=gassian, data=Data_mackerel_use_Ireland_select)
#     m0 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m1 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m2 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m3 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m4 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + Release_month + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m5 <-  gam(EW_move/1000 ~ offset(duration) + s(length, by= Release_year) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m6 <-  gam(EW_move/1000 ~ offset(duration) + s(length, by= Release_year) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m7 <-  gam(EW_move/1000 ~ offset(duration) + s(Longitude, Latitude, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m8 <-  gam(EW_move/1000 ~ offset(duration) + s(Longitude, Latitude, by= Release_year) + s(length) , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m9 <-  gam(EW_move/1000 ~ offset(duration) + s(Longitude, Latitude, by= Release_year) + length , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m10 <- gam(EW_move/1000 ~ offset(duration) + s(Longitude, Latitude, by= Release_year) + s(length, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
#
#     AICtab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
#     mbest <- m3
#     plot(mbest, all.terms=TRUE, residuals=TRUE, pages=1)
#     sim_select <- simulateResiduals(fittedModel = mbest, n = 1000, integerResponse = FALSE, plot=FALSE)
#     plot(sim_select, rank=T, quantreg = TRUE)
#     testResiduals(sim_select)
#     par(mfrow=c(3,2))
#     plotResiduals(Data_mackerel_use_Ireland_select$Length, sim_select$scaledResiduals, main = "length")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
#     plotResiduals(Data_mackerel_use_Ireland_select$Latitude, sim_select$scaledResiduals, main = "Latitude")
#     plotResiduals(Data_mackerel_use_Ireland_select$Longitude, sim_select$scaledResiduals, main = "Longitude")
#
#   # Movement North-South
#     m0NS <- gam(NS_move ~ s(length) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m1NS <- gam(NS_move ~ s(length) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m2NS <- gam(NS_move ~ s(length) + s(Longitude, Latitude) + s(duration), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m3NS <- gam(NS_move ~ s(length) + s(duration) + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m4NS <- gam(NS_move ~ s(length) + Release_month + s(Longitude, Latitude) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m5NS <- gam(NS_move ~ s(length, by= Release_year) + Release_month + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m6NS <- gam(NS_move ~ s(length, by= Release_year) + s(duration) + s(Longitude, Latitude), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m7NS <- gam(NS_move ~ s(Longitude, Latitude, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m8NS <- gam(NS_move ~ s(Longitude, Latitude, by= Release_year) + s(length) , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m9NS <- gam(NS_move ~ s(Longitude, Latitude, by= Release_year) + length , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m10NS <- gam(NS_move ~ s(Longitude, Latitude, by= Release_year) + s(length, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
#
#     AICtab(m0NS,m1NS,m2NS,m3NS,m4NS,m5NS,m6NS,m7NS,m8NS,m9NS,m10NS)
#     plot(m3NS, all.terms=TRUE, residuals=TRUE, pages=1)
#     sim_select <- simulateResiduals(fittedModel = m8, n = 1000, integerResponse = FALSE, plot=FALSE)
#     plot(sim_select, rank=T, quantreg = TRUE)
#     testResiduals(sim_select)
#     par(mfrow=c(3,2))
#     plotResiduals(Data_mackerel_use_Ireland_select$Length, sim_select$scaledResiduals, main = "length")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
#     plotResiduals(Data_mackerel_use_Ireland_select$Latitude, sim_select$scaledResiduals, main = "Latitude")
#     plotResiduals(Data_mackerel_use_Ireland_select$Longitude, sim_select$scaledResiduals, main = "Longitude")
#
#
# ##


