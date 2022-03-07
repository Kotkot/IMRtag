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

### Mvt rate = response variable in the analysis below
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
North_Ireland <- unique(Data_mackerel_all$ICES_Rectangle)[unlist(sapply(37:52, function(x) grep(x, unique(Data_mackerel_all$ICES_Rectangle))))]
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

### Focus on release around Ireland (not Iceland nor Bergen)
Data_mackerel_use_Ireland <- subset(Data_mackerel_final, Tag_area %in% c("South_Ireland", "North_Ireland"))
### Focus on release around Iceland (not Iceland nor Bergen)
Data_mackerel_use_Iceland <- subset(Data_mackerel_final, Tag_area %in% c("Iceland"))

### Do some explanatory analysis to explore what is causing these patterns of residuals
mean_tag_area <- Data_mackerel_use_Ireland_select %>% group_by(Tag_area) %>% summarize(median = median(log_rate))
mean_diff_tag_area <- as.numeric(abs(mean_tag_area[1,2] - mean_tag_area[2,2]))


########## Preparing the final data from Ireland for the analysis

Data_mackerel_use_Ireland1 <- Data_mackerel_use_Ireland

Adding_var <- function(data = Data_mackerel_use_Ireland) {
  data$Catch_month <- as.numeric(as.character(data$Catch_month))
  data$Catch_year <- as.numeric(as.character(data$Catch_year))
  data[which(data$Catch_month %in% c(1,2)),'Catch_month'] <- data[which(data$Catch_month %in% c(1,2)),'Catch_month'] + 12
  data[which(data$Catch_month %in% c(13,14)),'Catch_year'] <- data[which(data$Catch_month %in% c(13,14)),'Catch_year'] - 1

  ## Deriving data-frame for the 3 different migration cycles
  out1 <- data %>%
    filter(Catch_month %in% c(7:14), as.numeric(as.character(Catch_year)) == as.numeric(as.character(Release_year)), Release_year%in%2014:2020)
  out2 <- data %>%
    filter(Catch_month %in% c(7:14), as.numeric(as.character(Catch_year)) == as.numeric(as.character(Release_year))+1, Release_year%in%2014:2020)
  out3 <- data %>%
    filter(Catch_month %in% c(7:14),  as.numeric(as.character(Catch_year)) == as.numeric(as.character(Release_year))+2, Release_year%in%2014:2020)
  out1$Catch_year <- as.factor(out1$Catch_year)
  out2$Catch_year <- as.factor(out2$Catch_year)
  out3$Catch_year <- as.factor(out3$Catch_year)

  return(list(out1,out2,out3))

}


bla <- Adding_var(data = Data_mackerel_use_Ireland)
Data_mackerel_use_Ireland_select_origin <- bla[[1]]
Data_mackerel_use_Ireland_select_origin_year1 <- bla[[2]]
Data_mackerel_use_Ireland_select_origin_year2 <- bla[[3]]


bla <- Adding_var(data = Data_mackerel_use_Iceland)
Data_mackerel_use_Iceland_select_origin <- bla[[1]]
Data_mackerel_use_Iceland_select_origin_year1 <- bla[[2]]
Data_mackerel_use_Iceland_select_origin_year2 <- bla[[3]]

