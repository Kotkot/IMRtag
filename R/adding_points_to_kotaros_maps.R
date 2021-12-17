rm(list=ls())
library(tidyverse)
library(taggart)
library(lubridate)

Catch_data <- taggart::tg_catches(lower=T)
Mark_recap_data <- taggart::tg_expeditions(lower=T)
d <- Mark_recap_data %>% filter(year(releasedate) <2021)
d <- filter(d,
            ices_rectangle %in% unique(d$ices_rectangle)[unlist(sapply(31:52, function(x) grep(x, unique(d$ices_rectangle))))])

d <- filter(d,
            !(ices_rectangle %in% unique(d$ices_rectangle)[unlist(sapply("F", function(x) grep(x, unique(d$ices_rectangle))))]))


# Releases:
df_rel <- d %>% filter(between(year(releasedate),2014,2020))
df_rel <- mutate(df_rel,
                 iceslon =mapplots::ices.rect(ices_rectangle)$lon,
                 iceslat = mapplots::ices.rect(ices_rectangle)$lat)


# Recaptures:
recap_raw <- df_rel %>% filter(!is.na(recapturedate)) %>%
  left_join(Catch_data, by = "catchid") %>%
  mutate(RecaptureYear = ifelse(month(recapturedate)>2,
                                year(recapturedate),
                                year(recapturedate)-1),
         ReleaseYear = year(releasedate),
         YearsOut = RecaptureYear - ReleaseYear) %>%
  filter(YearsOut<3, between(ReleaseYear,2014,2020),
         !is.na(ices_rectangle.y))

# Get map:
xlim = range(c(recap_raw$clon,df_rel$longitude) , na.rm = T)
ylim = range(c(recap_raw$clat,df_rel$latitude), na.rm = T)
m <- ggplot2::map_data("world",
                       xlim = xlim + c(-.4,1),
                       ylim = ylim +c(-.1,.1) )

# Simple plot:
ggplot()+
  geom_polygon(data = m, aes(long, lat, group = group), fill = "grey") +
  coord_quickmap(xlim = xlim + c(-.5,.5), ylim = ylim +c(0,0)) +
  theme_bw()+
  # raw releases:
  geom_point(data = df_rel %>% distinct(longitude,latitude),    aes(x = longitude, y = latitude), col = "red", size = .1)+
  # ices rectangles:
  #geom_point(data = df_rel %>% distinct(iceslon,iceslat),    aes(x = iceslon, y = iceslat), col = "red", size = .1)+
  geom_point(data = recap_raw %>% distinct(clon,clat), aes(x = clon, y = clat), col = "#FEE724", size = .1)+
  theme(panel.background = element_rect(fill = "#00008B"),
        panel.grid = element_blank())
