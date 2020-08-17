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

### Installing required packages
	devtools::install_github("fishvice/taggart", dependencies = FALSE)

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

### Downloading data and checking it
	tg_catches()         %>% glimpse()
	tg_catches_bio()     %>% glimpse()
	tg_expeditions()     %>% glimpse()
	tg_expeditions_bio() %>% glimpse()

	Link_catch_sample <- jsonlite::fromJSON(paste0("http://smartfishsvc.hi.no/api/data/BioSamplesCatches/", species[1]))

	Catch_data <- tg_catches()
	Catch_bio <- tg_catches_bio()
	Mark_recap_data <- tg_expeditions()
	Mark_recap_bio <- tg_expeditions_bio()

	Data_mackerel <- Mark_recap_data %>% subset(., species == "mackerel")

	map_area <- get_map(location= c(left = -35, bottom = 51, right = 21, top = 72))
	#ggmap(map_area) + geom_point(data=Data_mackerel, aes(x=lo, y=la))

	view(Catch_data[which(Catch_data$reference_plant == "0f6a96b6-1b15-4fd2-ba2f-0d2aeae100bf"),])
	Catch_bio[which(Catch_bio$id == "e0d58869-c196-4fb5-a1c0-cba20f700da2"),]

### Define the data frame: creating and merging variables
	## Some base variables
		Data_mackerel_all <- Data_mackerel %>% left_join(Catch_data, by=c("catchid" = "pkid"))
		Data_mackerel_all$dist <- pointDistance(matrix(cbind(Data_mackerel_all$lo, Data_mackerel_all$la),ncol=2), matrix(cbind(Data_mackerel_all$cLon, Data_mackerel_all$cLat), ncol=2), lonlat=TRUE)
		Data_mackerel_all$duration <- Data_mackerel_all$recapturedate - Data_mackerel_all$relesedate
		Data_mackerel_all$Nb_years <- as.numeric(lubridate::year(as.Date(Data_mackerel_all$recapturedate))) - (lubridate::year(as.Date(Data_mackerel_all$relesedate))) + 1
		Data_mackerel_all$tagger <- as.factor(Data_mackerel_all$tagger )
		Data_mackerel_all$duration <- as.numeric(Data_mackerel_all$duration )
		Data_mackerel_all$duration_std <- Data_mackerel_all$duration/Data_mackerel_all$Nb_years
		Data_mackerel_all$duration_std1 <- Data_mackerel_all$duration%%365
		Data_mackerel_all$length <- as.numeric(Data_mackerel_all$length )
		Data_mackerel_all$la <- as.numeric(Data_mackerel_all$la )
		Data_mackerel_all$lo <- as.numeric(Data_mackerel_all$lo )
		Data_mackerel_all$dist <- as.numeric(Data_mackerel_all$dist )
		Data_mackerel_all$Release_year <- as.factor(lubridate::year(as.Date(Data_mackerel_all$relesedate)))
		Data_mackerel_all$Release_month <- as.factor(lubridate::month(as.Date(Data_mackerel_all$relesedate)))
		Data_mackerel_all$Release_month <- factor(Data_mackerel_all$Release_month, labels=c("05","06","08","09","10"))
		Data_mackerel_all$pos <- numFactor(Data_mackerel_all$lo, Data_mackerel_all$la)
		Data_mackerel_all$Catch_year <- as.factor(lubridate::year(as.Date(Data_mackerel_all$recapturedate)))
		Data_mackerel_all$Catch_month <- as.factor(lubridate::month(as.Date(Data_mackerel_all$recapturedate)))
		Data_mackerel_all$Release_day <- as.factor(lubridate::day(as.Date(Data_mackerel_all$relesedate)))
		Data_mackerel_all$Release_day <- factor(Data_mackerel_all$Release_day, labels=formatC(sort(unique(Data_mackerel_all$Release_day)), width=2, flag=0))
		Data_mackerel_all$Release_monthday <- paste(Data_mackerel_all$Release_month, Data_mackerel_all$Release_day, sep="_")
		Data_mackerel_all$Release_monthday <- ordered(Data_mackerel_all$Release_monthday)
		Data_mackerel_all$ID <- factor(rep(1, nrow(Data_mackerel_all)))
		Data_mackerel_all$firstyear <- as.Date(paste0(Data_mackerel_all$Release_year,"-01-01"))
		Data_mackerel_all$date <- as.numeric(as.Date(Data_mackerel_all$relesedate) - Data_mackerel_all$firstyear)

		## Mvt rate = response variable in the analysis below
		Data_mackerel_all$rate_annual <- Data_mackerel_all$dist / Data_mackerel_all$duration_std
		Data_mackerel_all$log_rate_annual <- log(Data_mackerel_all$dist/Data_mackerel_all$duration_std)
		Data_mackerel_all$rate_annual1 <- Data_mackerel_all$dist / Data_mackerel_all$duration_std1
		Data_mackerel_all$log_rate_annual1 <- log(Data_mackerel_all$dist/Data_mackerel_all$duration_std1)
		Data_mackerel_all$rate <- Data_mackerel_all$dist / Data_mackerel_all$duration
		Data_mackerel_all$log_rate <- log(Data_mackerel_all$dist/Data_mackerel_all$duration)

	## Mvt N-S or E-W
		Data_mackerel_all$EW_move <- as.numeric(pointDistance(matrix(cbind(Data_mackerel_all$lo, Data_mackerel_all$la),ncol=2), matrix(cbind(Data_mackerel_all$cLon, Data_mackerel_all$la), ncol=2), lonlat=TRUE))
		Data_mackerel_all$NS_move <- as.numeric(pointDistance(matrix(cbind(Data_mackerel_all$lo, Data_mackerel_all$la),ncol=2), matrix(cbind(Data_mackerel_all$lo, Data_mackerel_all$cLat), ncol=2), lonlat=TRUE))

	## Length bin
		Data_mackerel_all$length_bin <- cut(Data_mackerel_all$length, breaks=c(0, 33, 38, 45))
		Data_mackerel_all$length_bin  <- factor(Data_mackerel_all$length_bin, labels=c("(0,33cm]", "(33cm,38cm]", "(38cm,45cm]"))

	## Need to make some geographical division of the data
		West_Ireland <- unique(Data_mackerel_all$ices)[unlist(sapply(31:36, function(x) grep(x, unique(Data_mackerel_all$ices))))]
		North_Ireland <- unique(Data_mackerel_all$ices)[unlist(sapply(37:46, function(x) grep(x, unique(Data_mackerel_all$ices))))]
		Bergen <- unique(Data_mackerel_all$ices)[grep("F", unique(Data_mackerel_all$ices))]
		Iceland <- unique(Data_mackerel_all$ices)[grep(58, unique(Data_mackerel_all$ices))]
		Data_mackerel_all <- Data_mackerel_all %>% mutate(Tag_area = ifelse(ices %in% West_Ireland, "West_Ireland",
														 ifelse(ices %in% North_Ireland, "North_Ireland",
																ifelse(ices %in% Bergen, "Bergen",
																	   ifelse(ices %in% Iceland, "Iceland", "Weird")))))
		Data_mackerel_all <- Data_mackerel_all %>% mutate(Tag_area_large = ifelse(ices %in% c(West_Ireland,North_Ireland), "Ireland",
																					 ifelse(ices %in% Bergen, "Bergen",
																							ifelse(ices %in% Iceland, "Iceland", "Weird"))))

		Data_mackerel_all$Tag_area <- as.factor(Data_mackerel_all$Tag_area)
		Data_mackerel_all$Tag_area_large <- as.factor(Data_mackerel_all$Tag_area_large)

		Data_mackerel_all$Release_timing <- ifelse(Data_mackerel_all$Release_monthday < "05_15", "Early", "Late")
		Data_mackerel_all$Release_timing[which(Data_mackerel_all$Release_monthday >= "05_15" & Data_mackerel_all$Release_monthday <= "05_31")] <- "Mid"

		Data_mackerel_all$Category <- paste(Data_mackerel_all$Tag_area, Data_mackerel_all$Release_timing, sep="_")

### If I want to look at both the mark and recapture
	Data_mackerel_final <- Data_mackerel_all %>% subset(., !is.na(dist))

### Some visual inspection of the tagging data
	## How often do mackerel move south
		ggplot(Data_mackerel_all[which((Data_mackerel_all$la > Data_mackerel_all$cLat) ==TRUE),]) +
		geom_point(aes(x=lo, y=la)) + geom_point(aes(x=cLon, y=cLat), col="red")

	## Something is going on.... let's dig a bit more
		# Release month and location
			with(Data_mackerel_final, table(Release_month, Tag_area_large))
			with(Data_mackerel_final, table(Release_month, Tag_area_large, Release_year))
			with(Data_mackerel_final, table(Release_month, Tag_area))
			with(Data_mackerel_final, table(Release_month, Tag_area, Release_year))

		# Catch month and location
			with(Data_mackerel_final, table(Catch_month, ices_rectangle))
			with(Data_mackerel_final, table(Catch_month, recatch_ices_rectangle, Catch_year))
			with(Data_mackerel_final, table(Catch_month, Tag_area))
			with(Data_mackerel_final, table(Catch_month, Tag_area, Release_year))

		# Release month and length
			ggplot(Data_mackerel_final, aes(x=Release_month, y=length), alpha=0.5) + geom_boxplot() +
			facet_grid(Tag_area_large~ Release_year) + theme_bw() #+ coord_flip()

		# the duration
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

### Some plotting of the tagging location and recapture
### Maps + centre of gravity
	## By Year
		Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year) %>%
								mutate(mean_lon_tag=mean(lo),
									   mean_lat_tag=mean(la),
									   mean_lon_recap=mean(cLon),
									   mean_lat_recap=mean(cLat))
		Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year, length_bin) %>%
								mutate(mean_lon_recap_size=mean(cLon),
									   mean_lat_recap_size=mean(cLat))
		Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year, length_bin, Tag_area_large) %>%
								mutate(mean_lon_recap_size_area=mean(cLon),
									   mean_lat_recap_size_area=mean(cLat))
		Data_mackerel_final <- Data_mackerel_final %>% group_by(Release_year, Tag_area_large) %>%
								mutate(mean_lon_tag_area=mean(lo),
									   mean_lat_tag_area=mean(la))

		Data_mackerel_final_sameyear_sameyear <- Data_mackerel_final_sameyear[which(as.numeric(Data_mackerel_final_sameyear$Release_year) == as.numeric(Data_mackerel_final_sameyear$Catch_year)),]
		Data_mackerel_final_sameyear <- Data_mackerel_final_sameyear %>% group_by(Release_year) %>%
		  mutate(mean_lon_tag=mean(lo),
				 mean_lat_tag=mean(la),
				 mean_lon_recap=mean(cLon),
				 mean_lat_recap=mean(cLat))
		Data_mackerel_final_sameyear <- Data_mackerel_final_sameyear %>% group_by(Release_year, length_bin) %>%
		  mutate(mean_lon_recap_size=mean(cLon),
				 mean_lat_recap_size=mean(cLat))
		Data_mackerel_final_sameyear <- Data_mackerel_final_sameyear %>% group_by(Release_year, length_bin, Tag_area_large) %>%
		  mutate(mean_lon_recap_size_area=mean(cLon),
				 mean_lat_recap_size_area=mean(cLat))
		Data_mackerel_final_sameyear <- Data_mackerel_final_sameyear %>% group_by(Release_year, Tag_area_large) %>%
		  mutate(mean_lon_tag_area=mean(lo),
				 mean_lat_tag_area=mean(la))

		gg <- list()
		for (ijk in seq_along(2011:2019))
		{
		  gg[[ijk]] <- ggmap(map_area) + geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=lo, y=la)) +
			+ ggtitle((2011:2019)[ijk]) + theme(plot.title = element_text(hjust = 0.5)) +
			geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=cLon, y=cLat, col=length_bin), alpha=0.5) +
			geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_recap_size, y=mean_lat_recap_size, col=length_bin), shape=7, size=5) +
			geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=mean_lon_tag, y=mean_lat_tag), col="black", shape=7, size=5)

		}

		grid.arrange(gg[[1]],gg[[2]],gg[[3]],
					 gg[[4]],gg[[5]],gg[[6]],
					 gg[[7]],gg[[8]],gg[[9]],
					 nrow=3, ncol=3)

    ## by year, region and size
		# all years
        ggg <- list()
        for (ijk in seq_along(2011:2019))
        {
          ggg[[ijk]] <- ggmap(map_area) + geom_point(data=subset(Data_mackerel_final, Release_year==(2011:2019)[ijk]), aes(x=lo, y=la)) +
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

		# only 1 year subset
        gggg <- list()
        for (ijk in seq_along(2011:2019))
        {
          gggg[[ijk]] <- ggmap(map_area) + geom_point(data=subset(Data_mackerel_final_sameyear, Release_year==(2011:2019)[ijk]), aes(x=lo, y=la)) +
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


### Calculate linear distance between mark and recapture location
### Then standardize by number of days in between
### Model this with respect to the different covariates

	## Analysis over the whole area: but area west&north of Ireland is the only place with consistent release over the years

		## Using duration (not standardized)
		  # m01 <- gam(log_rate ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_final)
		  # m02 <- gam(log_rate ~ s(length) + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_final)
		  # m03 <- gam(log_rate ~ s(length) + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_final)
		  # m04 <- gam(rate ~ s(length) + s(lo, la) + Release_year, family=gaussian(link="log"), data=Data_mackerel_final)
		  # m05 <- gam(rate ~ s(length) + s(lo, la) + Release_year, family=Gamma(link = "log"), data=Data_mackerel_final)
		  # AICtab(m01,m02,m03,m04,m05)
		  # plot(m03, all.terms=TRUE, residuals=TRUE, pages=1)
		  # summary(m03)

		## Using duration (not standardized) but focusing only on within year recapture
		  Data_mackerel_use_select <- subset(Data_mackerel_final, duration < 365)
		  Data_mackerel_use_select <- Data_mackerel_final_sameyear
		  m0 <- gam(log_rate ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_select)
		  m1 <- gam(log_rate ~ s(length) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_select)
		  m2 <- gam(log_rate ~ s(length) + Release_month + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_select)
		  m3 <- gam(log_rate ~ s(length) + Release_month + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_select)
		  m4 <- gam(log_rate ~ s(length, by= Release_year) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_select)
		  m5 <- gam(log_rate ~ s(lo, la, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_select)
		  m6 <- gam(log_rate ~ s(lo, la, by= Tag_area) + Release_year + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_select)
		  AICtab(m1,m2,m3,m4,m5,m6)
		  plot(m5, all.terms=TRUE, residuals=TRUE, pages=1)
		  sim_select <- simulateResiduals(fittedModel = m5, n = 1000, integerResponse = FALSE, plot=FALSE)
		  plot(sim_select, rank=T, quantreg = TRUE)
		  testResiduals(sim_select)
		  par(mfrow=c(3,2))
		  plotResiduals(Data_mackerel_use_select$length, sim_select$scaledResiduals, main = "length")
		  plotResiduals(Data_mackerel_use_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
		  plotResiduals(Data_mackerel_use_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
		  plotResiduals(Data_mackerel_use_select$la, sim_select$scaledResiduals, main = "Latitude")
		  plotResiduals(Data_mackerel_use_select$lo, sim_select$scaledResiduals, main = "Longitude")

		## Using duration_std
		  # mm1 <- gam(log_rate_annual ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_final)
		  # mm2 <- gam(log_rate_annual ~ s(length) + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_final)
		  # mm3 <- gam(log_rate_annual ~ s(length) + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_final)
		  # mm4 <- gam(rate_annual ~ s(length) + s(lo, la) + Release_year, family=gaussian(link="log"), data=Data_mackerel_final)
		  # mm5 <- gam(rate_annual ~ s(length) + s(lo, la) + Release_year, family=Gamma(link = "log"), data=Data_mackerel_final)
		  # mm6 <- gam(rate_annual ~ s(length) + s(lo, la, by=Tag_area) + Release_year, family=Gamma(link = "log"), data=Data_mackerel_final)
		  # mm0 <- gam(log_rate_annual ~ s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_final)
		  #
		  # AICtab(mm0, mm1,mm2,mm3,mm4,mm5)
		  # plot(mm3, all.terms=TRUE, residuals=TRUE, pages=1)

		## Using duration std1
		  # mmm0 <- gam(log_rate_annual1 ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_final)
		  # mmm1 <- gam(log_rate_annual1 ~ s(length) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_final)
		  # mmm2 <- gam(log_rate_annual1 ~ s(length) + Release_month + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_final)
		  # mmm3 <- gam(log_rate_annual1 ~ s(length) + Release_month + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_final)
		  # mmm4 <- gam(log_rate_annual1 ~ s(length, by= Release_year) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_final)
		  # mmm5 <- gam(log_rate_annual1 ~ s(lo, la, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_final)
		  # mmm6 <- gam(log_rate_annual1 ~ s(lo, la, by= Tag_area) + s(length) + Release_year + Release_month , family=gaussian, data=Data_mackerel_final)
		  # mmm7 <- gam(log_rate_annual1 ~ s(lo, la, by= Tag_area)  + s(length, by= Release_year) + Release_year + Release_month , family=gaussian, data=Data_mackerel_final)
		  # mmm7a <- gam(log_rate_annual1 ~ s(lo, la, by= Tag_area)  + s(length, by= Release_year) + Release_month , family=gaussian, data=Data_mackerel_final)
		  # mmm7b <- gam(log_rate_annual1 ~ s(lo, la, by= Tag_area)  + s(length, by= Release_year) + Release_year , family=gaussian, data=Data_mackerel_final)
		  # mmm7c <- gam(log_rate_annual1 ~ s(lo, la, by= Tag_area)  + s(length, by= Release_year) , family=gaussian, data=Data_mackerel_final)

		  # AICtab(mmm0, mmm1,mmm2,mmm3,mmm4,mmm5,mmm6,mmm7,mmm7a,mmm7b,mmm7c)
		  # plot(mmm7, all.terms=TRUE, residuals=TRUE)
		  # sim_std1 <- simulateResiduals(fittedModel = mmm4, n = 1000, integerResponse = FALSE, plot=FALSE)
		  # plot(sim_std1, rank=T, quantreg = TRUE)
		  # testResiduals(sim_std1)
		  # par(mfrow=c(3,2))
		  # plotResiduals(Data_mackerel_final$length, sim_std1$scaledResiduals, main = "length")
		  # plotResiduals(Data_mackerel_final$Release_month, sim_std1$scaledResiduals, main = "Release_month")
		  # plotResiduals(Data_mackerel_final$Release_year, sim_std1$scaledResiduals, main = "Release_year")
		  # plotResiduals(Data_mackerel_final$la, sim_std1$scaledResiduals, main = "Latitude")
		  # plotResiduals(Data_mackerel_final$lo, sim_std1$scaledResiduals, main = "Longitude")

		#m1a <- gamm(rate ~ s(length) + s(lo, la), random=list(tagger= ~1), family=Gamma, data=Data_mackerel_final)
		#m2a <- gamm(rate ~ s(length) + s(lo, la) + s(date), random=list(tagger= ~1), family=Gamma, data=Data_mackerel_final)
		#mm1 <- glmmTMB(log_rate ~ length + mat(pos+0|ID) + Release_year, family=gaussian, data=Data_mackerel_final)


	## Limiting the analysis to releases around Ireland:

      Data_mackerel_use_Ireland <- subset(Data_mackerel_final, Tag_area %in% c("West_Ireland", "North_Ireland"))

		# Using duration (not standardized) but focusing on within year recapture
			Data_mackerel_use_Ireland_select <- subset(Data_mackerel_use_Ireland, duration < 180 & Release_year%in%2014:2019)
			Data_mackerel_use_Ireland_select$Release_year <- as.factor(as.character(Data_mackerel_use_Ireland_select$Release_year))

			m01 <- gam(log_rate ~ s(length) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
			m02 <- gam(log_rate ~ s(length) + Category, family=gaussian, data=Data_mackerel_use_Ireland_select)
			m03 <- gam(log_rate ~ s(length) + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m04 <- gam(log_rate ~ s(length) + s(lo, la) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m05 <- gam(log_rate ~ length + s(lo, la, by= Release_year) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m06 <- gam(log_rate ~ length + s(lo, la) + Release_year + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)

			m0 <- gam(log_rate ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m1 <- gam(log_rate ~ s(length) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m2 <- gam(log_rate ~ s(length) + Release_month + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m3 <- gam(log_rate ~ s(length) + Release_month + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
			m4 <- gam(log_rate ~ s(length, by= Release_year) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m5 <- gam(log_rate ~ s(lo, la, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland_select)
			m6 <- gam(log_rate ~ s(lo, la, by= Release_year) + s(length) , family=gaussian, data=Data_mackerel_use_Ireland_select)
			m7 <- gam(log_rate ~ s(lo, la, by= Release_year) + length , family=gaussian, data=Data_mackerel_use_Ireland_select)
			m8 <- gam(log_rate ~ s(lo, la, by= Release_year) + s(length, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)

			m10 <- gam(log_rate ~ length + s(lo, la) + s(date, by=Release_year) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
			m11 <- gam(log_rate ~ length + s(lo, la) + s(date, by=Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
			m12 <- gam(log_rate ~ length + s(lo, la, by=Release_year) + s(date, by=Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)

			AICtab(m1,m2,m3,m4,m5,m6,m7,m8,
				   m0,m01,m02,m03,m04,m05,m10,m11,m12)
			bestm <- m12
			plot(bestm, all.terms=TRUE, residuals=TRUE, pages=1)
			sim_select <- simulateResiduals(fittedModel = bestm, n = 1000, integerResponse = FALSE, plot=FALSE)
			plot(sim_select, rank=T, quantreg = TRUE)
			testResiduals(sim_select)
			par(mfrow=c(3,2))
			plotResiduals(Data_mackerel_use_Ireland_select$length, sim_select$scaledResiduals, main = "length")
			plotResiduals(Data_mackerel_use_Ireland_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
			plotResiduals(Data_mackerel_use_Ireland_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
			plotResiduals(Data_mackerel_use_Ireland_select$la, sim_select$scaledResiduals, main = "Latitude")
			plotResiduals(Data_mackerel_use_Ireland_select$lo, sim_select$scaledResiduals, main = "Longitude")

		# Using duration std1
			# mmm0_Ireland <- gam(log_rate_annual1 ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm01_Ireland <- gam(log_rate_annual1 ~ s(length) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
			# mmm02_Ireland <- gam(log_rate_annual1 ~ s(length) + Category, family=gaussian, data=Data_mackerel_use_Ireland_select)
			# mmm03_Ireland <- gam(log_rate_annual1 ~ s(length) + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_Ireland_select)
			# mmm04_Ireland <- gam(log_rate_annual1 ~ s(length) + s(lo, la) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
			# mmm05_Ireland <- gam(log_rate_annual1 ~ length + s(lo, la, by= Release_year) + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
			# mmm06_Ireland <- gam(log_rate_annual1 ~ length + s(lo, la) + Release_year + s(date, by=Tag_area), family=gaussian, data=Data_mackerel_use_Ireland_select)
			# AICtab(mmm0_Ireland, mmm01_Ireland,mmm02_Ireland,mmm03_Ireland,mmm04_Ireland,mmm05_Ireland,mmm06_Ireland)

			# mmm1_Ireland <- gam(log_rate_annual1 ~ s(length) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm2_Ireland <- gam(log_rate_annual1 ~ s(length) + Release_month + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm3_Ireland <- gam(log_rate_annual1 ~ s(length) + Release_month + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm4_Ireland <- gam(log_rate_annual1 ~ s(length, by= Release_year) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm5_Ireland <- gam(log_rate_annual1 ~ s(lo, la, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm6_Ireland  <- gam(log_rate_annual1 ~ s(length) + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland)
			# mmm7_Ireland  <- gam(log_rate_annual1 ~ s(length, by= Release_year) + s(lo, la, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland)

			# AICtab(mmm0_Ireland, mmm1_Ireland,mmm2_Ireland,mmm3_Ireland,mmm4_Ireland,mmm5_Ireland,mmm6_Ireland,mmm7_Ireland)

			# best_m <- mmm05_Ireland
			# plot(best_m, all.terms=TRUE, residuals=TRUE, pages=1)
			# sim_std1_Ireland <- simulateResiduals(fittedModel = best_m, n = 1000, integerResponse = FALSE, plot=FALSE)
			# plot(sim_std1_Ireland, rank=T, quantreg = TRUE)
			# testResiduals(sim_std1_Ireland)
			# par(mfrow=c(3,2))
			# plotResiduals(Data_mackerel_use_Ireland$length, sim_std1_Ireland$scaledResiduals, main = "length")
			# plotResiduals(Data_mackerel_use_Ireland$Release_month, sim_std1_Ireland$scaledResiduals, main = "Release_month")
			# plotResiduals(Data_mackerel_use_Ireland$Release_year, sim_std1_Ireland$scaledResiduals, main = "Release_year")
			# plotResiduals(Data_mackerel_use_Ireland$la, sim_std1_Ireland$scaledResiduals, main = "Latitude")
			# plotResiduals(Data_mackerel_use_Ireland$lo, sim_std1_Ireland$scaledResiduals, main = "Longitude")

		# Conclusion:
		# Residual pattern does not look good

        # Do some explanatory analysis to explore what is causing these patterns of residuals
			with(Data_mackerel_use_Ireland_select, plot(length, log_rate))
			with(subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland"), plot(length, log_rate))
			with(subset(Data_mackerel_use_Ireland_select, Tag_area=="West_Ireland"), plot(length, log_rate))
			with(subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland"), plot(date, log_rate))
			with(subset(Data_mackerel_use_Ireland_select, Tag_area=="West_Ireland"), plot(date, log_rate))
			ggplot(Data_mackerel_use_Ireland_select, aes(x=date, y=log_rate, col=Release_year)) + geom_point() + facet_grid(~Tag_area)
			mean_tag_area <- Data_mackerel_use_Ireland_select %>% group_by(Tag_area) %>% summarize(median = median(log_rate))
			mean_diff_tag_area <- as.numeric(abs(mean_tag_area[1,2] - mean_tag_area[2,2]))
			ggplot(Data_mackerel_use_Ireland_select, aes(x=date, y=log_rate, col=Release_year)) + geom_point() + facet_grid(~Tag_area) + geom_hline(yintercept = mean_diff_tag_area+ 8.7)+ geom_hline(yintercept =  8.7)

			brr <- subset(Data_mackerel_use_Ireland_select, Tag_area=="West_Ireland")
			brr$col <- ifelse(brr$log_rate<9.1,"red", ifelse(brr$log_rate>9.7, "blue", "purple"))
			ggmap(map_area) + geom_point(data=brr, aes(x=cLon, y=cLat, col=col), size=1)
			brr1 <- subset(Data_mackerel_use_Ireland_select, Tag_area=="North_Ireland")
			brr1$col <- ifelse(brr1$log_rate<8.7,"red", "purple")
			ggmap(map_area) + geom_point(data=brr1, aes(x=cLon, y=cLat, col=col), size=1)

			m0 <- gam(log_rate ~ s(length) + s(lo, la), family=gaussian, data=subset(brr, col=="red"))
			simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)
			m0 <- gam(log_rate ~ s(length) + s(lo, la), family=gaussian, data=subset(brr, col=="blue"))
			simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)
			m0 <- gam(log_rate ~ s(length) + s(lo, la), family=gaussian, data=subset(brr, col=="purple"))
			simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)
			m0 <- gam(log_rate ~ length + s(lo, la), family=gaussian, data=subset(brr1, col=="purple"))
			simulateResiduals(fittedModel = m0, n = 1000, integerResponse = FALSE, plot=TRUE)

        # This makes me think that I might need to develop a changepoint model (with K components)
		# I sometimes call it "mixture" model below but it is not a mixture model but a changepoint model
		# A Bayesian change point model has therefore been developped below
			test <- Data_mackerel_use_Ireland_select
			test <- test[order(test$log_rate),]

			m1 <- lm(log_rate ~ factor(Tag_area) + factor(Release_timing) + length:factor(Release_year), data=Data_mackerel_use_Ireland_select)
			m_frame <- model.frame(m1)
			XX <- model.matrix(m1, m_frame)

			# these are the change points values (not using the whole data to speed up the computing time)
			yyy <- Data_mackerel_use_Ireland_select$log_rate
			threshold_vals <-  as.numeric(quantile(yyy, seq(0.1, 0.9, by=0.05)))
			threshold_vals_group <- cut(test$log_rate, c(0, threshold_vals, 100))
			threshold_vals_group <- as.numeric(as.character(factor(threshold_vals_group, labels=1:18) ))
			threshold_vals_group_start <- c(1, which(diff(threshold_vals_group, lag=1) == 1)+1)
			threshold_vals_group_end <- c(which(diff(threshold_vals_group, lag=1) == 1), length(threshold_vals_group))

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
				  is_from_west=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "West_Ireland",1,0),
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
				  is_from_west=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "West_Ireland",1,0),
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


		## Now doing the same model but using TMB
				library(TMB)
				library(TMBhelper)
				use_version <- paste0(getwd(), "/src/mackerel_mvt_model")
				compile(paste0(use_version, ".cpp"))
				dyn.load(use_version)

				N_threshold <- 2
				data_tmb <- list(K=N_threshold,  # number of mixture components
            				     N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
            				     X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
            				     Nthres=length(threshold_vals),
            				     thresh=threshold_vals,
            				     mean_diff_tag_area= mean_diff_tag_area,
            				     is_from_west=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "West_Ireland",1,0),
            				     y = Data_mackerel_use_Ireland_select$log_rate
            				     )

				parameters_tmb <- list(beta = matrix(c(rep(10,3),runif(9,-2,2),runif(6*3,-0.05,0.05)),byrow=T, ncol=3),
				                       log_sigma = rep(log(0.2),N_threshold)
				                       )

				Map = list()

				op <- getwd()
				setwd(paste0(getwd(),"/src"))
				obj <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
				setwd(op)

				set.seed(1)
				rm(opt)
				opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
				sd_report_1break <- sdreport(obj1break)
				Check_Identifiable(opt1break)
				sigma <- as.vector(exp(summary(sd_report_1break, "fixed")[grep("log_sigma", rownames(summary(sd_report_1break, "fixed"))),1]))
				mu_pred <- matrix(summary(sd_report_1break, "report")[,1], ncol=2, byrow=FALSE) 
				
				# Calculating the actual prediction
				Likelihood <- matrix(NA, nrow=data_tmb$N, ncol=data_tmb$Nthres)
				for (n in 1:data_tmb$N){
					for (thr in 1:data_tmb$Nthres){
						if (data_tmb$y[n] < (data_tmb$thresh[thr]+data_tmb$is_from_west[n]*data_tmb$mean_diff_tag_area)){
							Likelihood[n,thr] = dnorm(data_tmb$y[n], mu_pred[n,1], sigma[1])
						}
						if (data_tmb$y[n] >= (data_tmb$thresh[thr]+data_tmb$is_from_west[n]*data_tmb$mean_diff_tag_area)){
							Likelihood[n,thr] = dnorm(data_tmb$y[n], mu_pred[n,2], sigma[2])
						}
					}
				}
					
				# The weighting factor is the average likelihood across the N datapoint per threshold value 
				weight <- apply(Likelihood, 2, mean)
				weight <- weight/sum(weight)
				plot(data_tmb$thresh, weight)
				
				Prediction <- rep(0, data_tmb$N)
				for (n in 1:data_tmb$N){
					for (thr in 1:data_tmb$Nthres){
						if (data_tmb$y[n] < (data_tmb$thresh[thr]+data_tmb$is_from_west[n]*data_tmb$mean_diff_tag_area)){
							Prediction[n] = Prediction[n]+weight[thr]*mu_pred[n,1]
						}
						if (data_tmb$y[n] >= (data_tmb$thresh[thr]+data_tmb$is_from_west[n]*data_tmb$mean_diff_tag_area)){
							Prediction[n] = Prediction[n]+weight[thr]*mu_pred[n,2]
						}
					}
				}
				
				plot(Prediction, data_tmb$y); abline(0,1)
				qqnorm(y=(Prediction-data_tmb$y)/sd(Prediction-data_tmb$y))
				abline(0,1, lty=2)
				
						

############ Some nice plots to lok at mark recapture pattern

  ## some quick investigation
    Data_mackerel_all %>% subset(Release_year==2018 &
                        Tag_area %in% c("West_Ireland", "North_Ireland")) %>%
                        group_by(Tag_area) %>% dplyr::select(Tag_area, Release_monthday) %>% table()

  ## Plot to investigate the location of recapture from Ireland, by year
    plot_Yr <- function(Year)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("West_Ireland", "North_Ireland"))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year &
                                                    Tag_area %in% c("West_Ireland", "North_Ireland"))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))

      ncount_release <- dat_release %>% group_by(Tag_area) %>% count()
      ncount_release$lon = c(-10, -14)
      ncount_release$lat = c(57, 51)
      ncount_release$label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(Tag_area) %>% count()
      ncount_recap$lon = c(-8, -8)
      ncount_recap$lat = c(69, 70)
      ncount_recap$label = paste0("n=", ncount_recap$n)

      hull <- dat_recap %>% dplyr::select(Tag_area, cLon, cLat) %>% group_by(Tag_area) %>%  dplyr::slice(chull(cLon, cLat))
      gravity <- hull %>% group_by(Tag_area) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Tag_area, mean_lon, mean_lat) %>% dplyr::distinct()


      g1 <- ggmap(map_area) + geom_point(data=subset(dat_release, Tag_area=="From W. Ireland"), aes(x=lo, y=la), size=0.3, col="red") +
        geom_point(data=subset(dat_recap, Tag_area=="From W. Ireland"), aes(x=cLon, y=cLat), size=0.3, col="red") +
        geom_text(data=subset(ncount_recap, Tag_area=="From W. Ireland"), aes(x=lon, y=lat, label=label), col="red") +
        geom_text(data=subset(ncount_release, Tag_area=="From W. Ireland"), aes(x=lon, y=lat, label=label), col="red")  +
        ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) +
        geom_point(data=subset(dat_release, Tag_area=="From N. Ireland"), aes(x=lo, y=la), size=0.3, col="blue") +
        geom_point(data=subset(dat_recap, Tag_area=="From N. Ireland"), aes(x=cLon, y=cLat), size=0.3, col="blue") +
        geom_text(data=subset(ncount_recap, Tag_area=="From N. Ireland"), aes(x=lon, y=lat, label=label), col="blue") +
        geom_text(data=subset(ncount_release, Tag_area=="From N. Ireland"), aes(x=lon, y=lat, label=label), col="blue") +
        geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = factor(Tag_area))) +
        scale_fill_manual(name = "Convex hull recapture", values=c("lightblue", "pink")) +
        geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        geom_point(data=gravity, aes(x=mean_lon, y=mean_lat, col=Tag_area), shape=13, size=4) +
        scale_color_manual(name = "Center of gravity", values=c("darkblue", "darkred"))

      g1
      return(g1)
    }

    ## 2011-2013, there is no recapture originating from both West and North side of Ireland
    p2014 <- plot_Yr(2014)
    p2015 <- plot_Yr(2015)
    p2016 <- plot_Yr(2016)
    p2017 <- plot_Yr(2017)
    p2018 <- plot_Yr(2018)
    p2019 <- plot_Yr(2019)

    ggg <- arrangeGrob(p2014,p2015,p2016,p2017,p2018,p2019, nrow=3, ncol=2)
    ggsave(ggg, filename=paste0("../plots/mark_recap_2014-2019.pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

  #### Plot to investigate the location of recapture from Ireland, by year and
  #### Time of release
  #### To Determine if what matter is time of release or duration
    plot_Yr_timerelease <- function(Year)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("West_Ireland", "North_Ireland"))
      dat_release$julian <-  as.numeric(julian(dat_release$relesedate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year &
                                                    Tag_area %in% c("West_Ireland", "North_Ireland"))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))
      # dat_recap <- dat_recap %>% mutate(Release_timing = ifelse(Release_monthday <= "5_14", "Early",
                                                                # ifelse((Release_monthday <= "5_31" & Release_monthday >= "5_15"), "Mid", "Late")))

      dat_recap$Release_timing <- ifelse(dat_recap$Release_monthday < "05_15", "Early", "Late")
      dat_recap$Release_timing[which(dat_recap$Release_monthday >= "05_15" & dat_recap$Release_monthday <= "05_31")] <- "Mid"

      ncount_release <- dat_release %>% group_by(Tag_area) %>% count()
      ncount_release$lon = c(-17, -17)
      ncount_release$lat = c(57, 51)
      ncount_release$label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(Tag_area) %>% count()
      ncount_recap$lon = c(-8, -8)
      ncount_recap$lat = c(69, 70)
      ncount_recap$label = paste0("n=", ncount_recap$n)

      dat_recap$Category <- paste(dat_recap$Tag_area, dat_recap$Release_timing, sep="_")

      Select_catch <- subset(Catch_data, subset=c(catchdate ==Year))
      Select_catch <- subset(Select_catch, subset=c(!is.na(cLon) | !is.na(cLat)))
      Select_catch <- subset(Select_catch, cLat<72)
      hull_catch <- Select_catch %>% dplyr::slice(chull(cLon, cLat)) %>% dplyr::select(cLon, cLat)

      hull <- dat_recap %>% dplyr::select(Tag_area, Release_timing, cLon, cLat) %>% group_by(Tag_area, Release_timing) %>%  dplyr::slice(chull(cLon, cLat))
      hull$Category <- paste(hull$Tag_area, hull$Release_timing, sep="_")
      hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
      #hull$Category <- as.factor(hull$Category)
      #hull$Category <- factor(hull$Category, levels=c("From N. Ireland_Early", "From N. Ireland_Mid",
      #                                                "From N. Ireland_Late", "From W. Ireland_Early",
      #                                                "From W. Ireland_Mid", "From W. Ireland_Late"))
      cats <- c("From N. Ireland_Early", "From N. Ireland_Mid",
                "From N. Ireland_Late", "From W. Ireland_Early",
                "From W. Ireland_Mid", "From W. Ireland_Late")
      gravity <- hull %>% group_by(Category) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
      hull_new <- as.data.frame(hull)
      coordinates(hull_new) <- c("cLon", "cLat")
      crs(hull_new) <- "+proj=longlat +datum=WGS84"
      new_proj <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
      hull_new <- spTransform(hull_new, new_proj)
      hull_new <- as.data.frame(hull_new)
      area <- sapply(unique(hull_new$Category), function(x) Polygon(subset(hull_new, Category==x)[,c('cLon','cLat')])@area)
      area_category <- rep(NA, length(cats))
      names(area_category) <- cats
      area_category[match(names(area), names(area_category))] <- area

      prop_size <- dat_recap %>% group_by(Category,length_bin) %>% summarize(n=n()) %>% mutate(prop = prop.table(n)) %>% filter(n>1)
      prop_size$r <- 1
      Ngroups <- length(unique(prop_size$Category))
      prop_size$lon <- -33
      prop_size$lat <- as.factor(prop_size$Category)
      prop_size$lat <- as.numeric(as.character(factor(prop_size$lat, labels=seq(58,56+2.5*length(levels(prop_size$lat)),by=2.5))))
      prop_size_wide <- spread(prop_size[,-3], length_bin, prop)
      prop_size_wide[is.na(prop_size_wide)] <- 0

      prop_size_text <- data.frame(lon=-31, lat=seq(58,56+2.5*length(unique(prop_size$lat)),by=2.5), label=unique(prop_size$Category))

      g1 <- ggmap(map_area) + geom_jitter(data=subset(dat_release, Tag_area=="From W. Ireland"), aes(x=lo, y=la,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_point(data=subset(dat_recap, Tag_area=="From W. Ireland"), aes(x=cLon, y=cLat), size=0.3, col="black") +
        geom_jitter(data=subset(dat_release, Tag_area=="From N. Ireland"), aes(x=lo, y=la,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_text(data=subset(ncount_recap, Tag_area=="From W. Ireland"), aes(x=lon, y=lat, label=label), col="red") +
        geom_text(data=subset(ncount_release, Tag_area=="From W. Ireland"), aes(x=lon, y=lat, label=label), col="red")  +
        geom_text(data=subset(ncount_recap, Tag_area=="From N. Ireland"), aes(x=lon, y=lat, label=label), col="blue") +
        geom_text(data=subset(ncount_release, Tag_area=="From N. Ireland"), aes(x=lon, y=lat, label=label), col="blue") +
        geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat, fill = Category)) +
        scale_color_viridis_c() +
        geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        geom_text(data=data.frame(x=-25,y=54,label="Nsamp release:"), aes(x=x, y=y, label=label), col="black", fontface="bold") +
        ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5)) +
        geom_polygon(data = hull_catch, alpha = 0.5, aes(x=cLon, y=cLat), col="black", fill=NA) +
        geom_scatterpie(aes(x=lon, y=lat, group = Category, r=r),
                        data = prop_size_wide, cols = colnames(prop_size_wide[,-c(1:4)]))+
        geom_text(data=prop_size_text, aes(x=lon, y=lat, label=label), size=2, hjust = 0)
      g1


      return(list(g1, area_category))
    }

    ## 2011-2013, there is no recapture originating from both West and North side of Ireland
    p2014_timerelease <- plot_Yr_timerelease(2014)
    p2015_timerelease <- plot_Yr_timerelease(2015)
    p2016_timerelease <- plot_Yr_timerelease(2016)
    p2017_timerelease <- plot_Yr_timerelease(2017)
    p2018_timerelease <- plot_Yr_timerelease(2018)
    p2019_timerelease <- plot_Yr_timerelease(2019)

    areas <- rbind(p2014_timerelease[[2]],
                   p2015_timerelease[[2]],
                   p2016_timerelease[[2]],
                   p2017_timerelease[[2]],
                   p2018_timerelease[[2]],
                   p2019_timerelease[[2]])
    areas.std <- apply(areas, 1, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
    colnames(areas.std) <- 2014:2019
    # Testing mid W. Ireland bigger than others
    t.test(areas[-4,5], areas[-4,1], alternative = "greater")
    t.test(areas[-4,5], areas[-4,2], alternative = "greater")
    t.test(areas[-4,5], areas[-4,3], alternative = "greater")
    t.test(areas[-4,5], areas[-4,4], alternative = "greater")
    t.test(areas[-4,5], areas[-4,6], alternative = "greater")
    # Testing mid N. Ireland smaller than others
    t.test(areas[,2], areas[,1], alternative = "less")
    t.test(areas[,2], areas[,3], alternative = "less")
    t.test(areas[,2], areas[,4], alternative = "less")
    t.test(areas[,2], areas[,5], alternative = "less")
    t.test(areas[,2], areas[,6], alternative = "less")

    t.test(areas[4,], areas[1,], alternative = "two.sided")
    t.test(areas[4,], areas[2,], alternative = "two.sided")
    t.test(areas[4,], areas[3,], alternative = "two.sided")
    t.test(areas[4,], areas[4,], alternative = "two.sided")
    t.test(areas[4,], areas[5,], alternative = "two.sided")

    #ggplot() + gg_text()
    ggg_timerelease <- arrangeGrob(p2014_timerelease[[1]],p2015_timerelease[[1]],p2016_timerelease[[1]],
                       p2017_timerelease[[1]],p2018_timerelease[[1]],p2019_timerelease[[1]],
                       nrow=3, ncol=2)
    ggsave(ggg_timerelease, filename=paste0("../plots/mark_recap_release_2014-2019.pdf"), width=32, height=32, units="cm", device = "pdf", dpi = 400)

  ### Plot to investigate the location of recapture from Ireland, by year and
  ### Duration between release and catch
  ### To Determine if what matter is time of release or duration
    plot_Yr_duration <- function(Year)
    {
      dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                                    Tag_area %in% c("West_Ireland", "North_Ireland"))
      dat_release$julian <-  as.numeric(julian(dat_release$relesedate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
      dat_release$Tag_area <- as.character(dat_release$Tag_area)
      dat_release$Tag_area <- as.factor(dat_release$Tag_area)
      dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))
      dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year &
                                                    Tag_area %in% c("West_Ireland", "North_Ireland"))
      dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
      dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)
      dat_recap$Tag_area <- factor(dat_recap$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))

      dat_recap$duration_bin <- cut(dat_recap$duration, breaks=c(0,100,150,250))

      ncount_release <- dat_release %>% group_by(Tag_area) %>% count()
      ncount_release$lon = c(-17, -17)
      ncount_release$lat = c(57, 51)
      ncount_release$label = paste0("n=", ncount_release$n)

      ncount_recap <- dat_recap %>% group_by(Tag_area) %>% count()
      ncount_recap$lon = c(-8, -8)
      ncount_recap$lat = c(69, 70)
      ncount_recap$label = paste0("n=", ncount_recap$n)


      hull <- dat_recap %>% dplyr::select(Tag_area, duration_bin, cLon, cLat) %>% group_by(Tag_area, duration_bin) %>%  dplyr::slice(chull(cLon, cLat))
      hull$Category <- paste(hull$Tag_area, hull$duration_bin, sep="_")
      hull <- hull %>% group_by(Category) %>% filter(n()> 2) %>% ungroup()
      gravity <- hull %>% group_by(Category) %>%
        mutate(mean_lon= mean(cLon), mean_lat= mean(cLat)) %>%
        dplyr::select(Category, mean_lon, mean_lat) %>% dplyr::distinct()
      #area <- hull %>% group_by(Category) %>% Polygon()

      g1 <- ggmap(map_area) + geom_jitter(data=subset(dat_release, Tag_area=="From W. Ireland"), aes(x=lo, y=la,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_point(data=subset(dat_recap, Tag_area=="From W. Ireland"), aes(x=cLon, y=cLat), size=0.3, col="black") +
        geom_jitter(data=subset(dat_release, Tag_area=="From N. Ireland"), aes(x=lo, y=la,col=julian), size=0.3, shape=4, width = 0.2, height = 0.2) +
        geom_text(data=subset(ncount_recap, Tag_area=="From W. Ireland"), aes(x=lon, y=lat, label=label), col="red") +
        geom_text(data=subset(ncount_release, Tag_area=="From W. Ireland"), aes(x=lon, y=lat, label=label), col="red")  +
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

    ## 2011-2013, there is no recapture originating from both West and North side of Ireland
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
      m0_Iceland <- gam(log_rate ~ length + s(lo, la), family=gaussian, data=Data_mackerel_use_Iceland_select)
      m1_Iceland <- gam(log_rate ~ length + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_Iceland_select)
      m2_Iceland <- gam(log_rate ~ length + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Iceland_select)
      m3_Iceland <- gam(log_rate ~ length + s(lo, la), family=gaussian, data=Data_mackerel_use_Iceland_select)
      m4_Iceland <- gam(log_rate ~ length + Release_year, family=gaussian, data=Data_mackerel_use_Iceland_select)
      m5_Iceland <- gam(log_rate ~ length, family=gaussian, data=Data_mackerel_use_Iceland_select)

      AICtab(m1_Iceland,m2_Iceland,m3_Iceland,m4_Iceland,m5_Iceland)
      plot(m4_Iceland, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_select_Iceland <- simulateResiduals(fittedModel = m4_Iceland, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_select_Iceland, rank=T, quantreg = TRUE)
      testResiduals(sim_select_Iceland)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Iceland_select$length, sim_select_Iceland$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Iceland_select$Release_month, sim_select_Iceland$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Iceland_select$Release_year, sim_select_Iceland$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Iceland_select$la, sim_select_Iceland$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Iceland_select$lo, sim_select_Iceland$scaledResiduals, main = "Longitude")

      ## Using duration std1
      mmm0_Iceland <- gam(log_rate_annual1 ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm1_Iceland <- gam(log_rate_annual1 ~ s(length) + s(lo, la) + s(date), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm2_Iceland <- gam(log_rate_annual1 ~ s(length) + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Iceland)
      mmm3_Iceland <- gam(log_rate_annual1 ~ s(length, by= Release_year) + s(lo, la), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm4_Iceland <- gam(log_rate_annual1 ~ s(lo, la, by= Release_year) + s(length), family=gaussian, data=Data_mackerel_use_Iceland)
      mmm5_Iceland <- gam(log_rate_annual1 ~ s(length) + Release_year, family=gaussian, data=Data_mackerel_use_Iceland)

      AICtab(mmm0_Iceland, mmm1_Iceland,mmm2_Iceland,mmm3_Iceland,mmm4_Iceland,mmm5_Iceland,mmm6_Iceland)
      plot(mmm3_Iceland, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_std1_Iceland <- simulateResiduals(fittedModel = mmm3_Iceland, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_std1_Iceland, rank=T, quantreg = TRUE)
      testResiduals(sim_std1_Iceland)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Iceland$length, sim_std1_Iceland$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Iceland$Release_month, sim_std1_Iceland$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Iceland$Release_year, sim_std1_Iceland$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Iceland$la, sim_std1_Iceland$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Iceland$lo, sim_std1_Iceland$scaledResiduals, main = "Longitude")



  #### Limiting the analysis to releases around Bergen:

      Data_mackerel_use_Bergen <- subset(Data_mackerel_final, Tag_area %in% c("Bergen"))

      ## Using duration std1
      mmm0_Bergen <- gam(log_rate_annual1 ~ 1, family=gaussian, data=Data_mackerel_use_Bergen)
      mmm1_Bergen <- gam(log_rate_annual1 ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_Bergen)
      mmm2_Bergen <- gam(log_rate_annual1 ~ s(lo, la), family=gaussian, data=Data_mackerel_use_Bergen)
      mmm3_Bergen <- gam(log_rate_annual1 ~ s(length), family=gaussian, data=Data_mackerel_use_Bergen)

      AICtab(mmm0_Bergen, mmm1_Bergen,mmm2_Bergen,mmm3_Bergen)
      plot(mmm1_Bergen, all.terms=TRUE, residuals=TRUE, pages=1)
      sim_std1_Bergen <- simulateResiduals(fittedModel = mmm1_Bergen, n = 1000, integerResponse = FALSE, plot=FALSE)
      plot(sim_std1_Bergen, rank=T, quantreg = TRUE)
      testResiduals(sim_std1_Bergen)
      par(mfrow=c(3,2))
      plotResiduals(Data_mackerel_use_Bergen$length, sim_std1_Bergen$scaledResiduals, main = "length")
      plotResiduals(Data_mackerel_use_Bergen$Release_month, sim_std1_Bergen$scaledResiduals, main = "Release_month")
      plotResiduals(Data_mackerel_use_Bergen$Release_year, sim_std1_Bergen$scaledResiduals, main = "Release_year")
      plotResiduals(Data_mackerel_use_Bergen$la, sim_std1_Bergen$scaledResiduals, main = "Latitude")
      plotResiduals(Data_mackerel_use_Bergen$lo, sim_std1_Bergen$scaledResiduals, main = "Longitude")



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
#   Data_mackerel_use_Ireland <- subset(Data_mackerel_final, Tag_area %in% c("West_Ireland", "North_Ireland"))
#
#   ## Using duration (not standardized) but focusing on within year recapture
#   Data_mackerel_use_Ireland_select <- subset(Data_mackerel_use_Ireland, duration < 200)
#   Data_mackerel_use_Ireland_select$location <- with(Data_mackerel_use_Ireland_select, glmmTMB::numFactor(lo,la))
#   levels(Data_mackerel_use_Ireland_select$location)
#   Data_mackerel_use_Ireland_select$group <- factor(rep(1, nrow(Data_mackerel_use_Ireland_select)))
#   # Movement East_West
#     library(glmmTMB)
#     #mm0 <- glmmTMB(EW_move/1000 ~ offset(duration) + length + Tag_area, family=gassian, data=Data_mackerel_use_Ireland_select)
#     m0 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + Tag_area, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m1 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m2 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m3 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m4 <-  gam(EW_move/1000 ~ offset(duration) + s(length) + Release_month + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m5 <-  gam(EW_move/1000 ~ offset(duration) + s(length, by= Release_year) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m6 <-  gam(EW_move/1000 ~ offset(duration) + s(length, by= Release_year) + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m7 <-  gam(EW_move/1000 ~ offset(duration) + s(lo, la, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m8 <-  gam(EW_move/1000 ~ offset(duration) + s(lo, la, by= Release_year) + s(length) , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m9 <-  gam(EW_move/1000 ~ offset(duration) + s(lo, la, by= Release_year) + length , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m10 <- gam(EW_move/1000 ~ offset(duration) + s(lo, la, by= Release_year) + s(length, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
#
#     AICtab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
#     mbest <- m3
#     plot(mbest, all.terms=TRUE, residuals=TRUE, pages=1)
#     sim_select <- simulateResiduals(fittedModel = mbest, n = 1000, integerResponse = FALSE, plot=FALSE)
#     plot(sim_select, rank=T, quantreg = TRUE)
#     testResiduals(sim_select)
#     par(mfrow=c(3,2))
#     plotResiduals(Data_mackerel_use_Ireland_select$length, sim_select$scaledResiduals, main = "length")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
#     plotResiduals(Data_mackerel_use_Ireland_select$la, sim_select$scaledResiduals, main = "Latitude")
#     plotResiduals(Data_mackerel_use_Ireland_select$lo, sim_select$scaledResiduals, main = "Longitude")
#
#   # Movement North-South
#     m0NS <- gam(NS_move ~ s(length) + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m1NS <- gam(NS_move ~ s(length) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m2NS <- gam(NS_move ~ s(length) + s(lo, la) + s(duration), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m3NS <- gam(NS_move ~ s(length) + s(duration) + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m4NS <- gam(NS_move ~ s(length) + Release_month + s(lo, la) + Release_year, family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m5NS <- gam(NS_move ~ s(length, by= Release_year) + Release_month + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m6NS <- gam(NS_move ~ s(length, by= Release_year) + s(duration) + s(lo, la), family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m7NS <- gam(NS_move ~ s(lo, la, by= Release_year) + s(length) + Release_month , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m8NS <- gam(NS_move ~ s(lo, la, by= Release_year) + s(length) , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m9NS <- gam(NS_move ~ s(lo, la, by= Release_year) + length , family=gaussian, data=Data_mackerel_use_Ireland_select)
#     m10NS <- gam(NS_move ~ s(lo, la, by= Release_year) + s(length, by= Release_year), family=gaussian, data=Data_mackerel_use_Ireland_select)
#
#     AICtab(m0NS,m1NS,m2NS,m3NS,m4NS,m5NS,m6NS,m7NS,m8NS,m9NS,m10NS)
#     plot(m3NS, all.terms=TRUE, residuals=TRUE, pages=1)
#     sim_select <- simulateResiduals(fittedModel = m8, n = 1000, integerResponse = FALSE, plot=FALSE)
#     plot(sim_select, rank=T, quantreg = TRUE)
#     testResiduals(sim_select)
#     par(mfrow=c(3,2))
#     plotResiduals(Data_mackerel_use_Ireland_select$length, sim_select$scaledResiduals, main = "length")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_month, sim_select$scaledResiduals, main = "Release_month")
#     plotResiduals(Data_mackerel_use_Ireland_select$Release_year, sim_select$scaledResiduals, main = "Release_year")
#     plotResiduals(Data_mackerel_use_Ireland_select$la, sim_select$scaledResiduals, main = "Latitude")
#     plotResiduals(Data_mackerel_use_Ireland_select$lo, sim_select$scaledResiduals, main = "Longitude")
#
#
# ##


