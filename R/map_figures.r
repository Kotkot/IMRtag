# rm(list =ls())
# library(tidyverse)
# library(taggart)
#
#
# Catch_data <- taggart::tg_catches()
# Catch_bio <- taggart::tg_catches_bio()
# Mark_recap_data <- taggart::tg_expeditions()
# Mark_recap_bio <- taggart::tg_expeditions_bio()
# ------------------------
library(sf)
ICESecoregion <- st_read("shapefile/ICES_ecoregions_20171207_erase_ESRI.shp")
ICESecoregion <- ICESecoregion %>% filter(Ecoregion %in% c("Icelandic Waters",
                                                           "Norwegian Sea",
                                                           "Greater North Sea",
                                                           "Celtic Seas",
                                                           "Faroes"
))
nc_sp <- as_Spatial(ICESecoregion$geom, IDs=ICESecoregion$Ecoregion)
#sp::plot(nc_sp)
countymap <- fortify(nc_sp) %>% filter(!hole)

Mark_recap_data$ReleaseYear <- lubridate::year(Mark_recap_data$ReleaseDate)
d <- Mark_recap_data %>% filter(ReleaseYear <=2021)
d <- filter(d,
            ICES_Rectangle %in% unique(d$ICES_Rectangle)[unlist(sapply(31:52, function(x) grep(x, unique(d$ICES_Rectangle))))])

d <- filter(d,
            !(ICES_Rectangle %in% unique(d$ICES_Rectangle)[unlist(sapply("F", function(x) grep(x, unique(d$ICES_Rectangle))))]))

dat_release_new <- as_tibble(d) %>%
  dplyr::select(ICES_Rectangle, ReleaseYear) %>% count(ICES_Rectangle,
                                                       ReleaseYear)
dat_release_new$lon <- mapplots::ices.rect(dat_release_new$ICES_Rectangle)$lon
dat_release_new$lat <- mapplots::ices.rect(dat_release_new$ICES_Rectangle)$lat
#dat_release_new <- filter(dat_release_new, Latitude > 40, Longitude < 20)

# Calculate range of longitudes and latitudes
xlim = range(dat_release_new$lon, na.rm = T)
ylim = range(dat_release_new$lat, na.rm = T)



# Get map:
m <- ggplot2::map_data("world",
                       xlim = xlim + c(-.4,1),
                       ylim = ylim +c(-.1,.1) )
dat_release_count <- group_by(dat_release_new, ReleaseYear) %>% summarize(N = sum(n))
# Create plot:
dat_release_new<- left_join(dat_release_new, dat_release_count, by = "ReleaseYear")
dat_release_new <- mutate(dat_release_new, YearN = paste0(ReleaseYear, ", N = ", N))
dat_release_new$n

dat_release_new <- dat_release_new %>% mutate(ncat = cut(n,
                                                         breaks = c(-Inf, 1000, 3000, 5000, Inf),
                                                         labels = c("<1k", "1-3k", "3-5k", ">5k")))

dat_release_new<-dat_release_new %>% filter(between(ReleaseYear,2014,2020))

# For adding density curve
recap_raw <- d %>% filter(!is.na(RecaptureDate)) %>% left_join(Catch_data, by = "CatchID") %>%
  mutate(RecaptureYear = ifelse(lubridate::month(RecaptureDate)>2, lubridate::year(RecaptureDate), lubridate::year(RecaptureDate)-1),
         YearsOut = RecaptureYear - ReleaseYear) %>% filter(YearsOut<3, between(ReleaseYear,2014,2020),
                                                            !is.na(ICES_Rectangle.y))

p1 <- list()
plength <- plength2 <- list()
for(i in 1:length(unique(dat_release_new$ReleaseYear))){
  RL = (2014:2020)[i]
  df <- filter(dat_release_new, ReleaseYear == RL)
  df.count <- summarize(df, N = sum(n)) %>% mutate(Nlab = paste0("N=",N), lon = -Inf, lat = Inf)
(p1[[i]] <- ggplot() +
  geom_polygon(data = m, aes(long, lat, group = group), fill = "grey") +
  geom_tile(data = df,
            aes(lon, lat, fill = ncat))+#,# col = CatchWeight)) +
  coord_quickmap(xlim = xlim + c(-.5,.5), ylim = ylim +c(0,0)) +
  facet_wrap( ~ReleaseYear, ncol = 1, strip.position = "left") +
  theme_bw()+
  #guides(fill = guide_colorbar(barheight = 35))+
  #scale_fill_gradientn(colours = heat.colors(50)[45:5],
 #                      name = "#Releases") +
  scale_fill_manual(values = c("yellow", "orange", "red", "brown3"),
                    name = "#Releases")+
  theme(strip.background = element_rect(fill = "white", color ="white"),
        text = element_text(size = 26),
        legend.position = "left",
        legend.direction = "horizontal",
        strip.placement = "outside",
        panel.grid = element_blank(), axis.title = element_blank(),
        plot.margin = unit(c(2,2,2,2), "mm"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA))+
  scale_x_continuous("Longitude", breaks = seq(-15,5,5),
                     labels = paste0(abs(seq(-15,5,5)), "°W"))+
  scale_y_continuous("Latitude", breaks = seq(50,70,2),
                     labels = paste0(seq(50,70,2), "°"))+
    geom_text(data = df.count, aes(label = Nlab, x = lon, y=lat), hjust = -0.05, vjust = 1.2,
              size = 9))
  plength[[i]] <- ggplot(data=filter(d, ReleaseYear == RL), aes(x = Length)) +
     # stat_bin(aes(y=..density..),fill = "orange", binwidth  =1) +
    geom_density(aes(col = "Releases",fill = "Releases"), bw = .5, lwd = 1.4, alpha = .1) +
#      geom_density(aes(x=Length,
 #                      y=..density..))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,55,5),
                       minor_breaks= 0:55,
                       limits = c(21.5,45.5)) +
    geom_vline(xintercept = 35, lty =2, col = "black")+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,.3,0.05), labels = paste0(seq(0,30,5), "%")) +
      theme_bw() +
    scale_color_manual(name = "", values = c("red", "orange"))+
    scale_fill_manual(name = "", values = c("red", "orange"))+
    geom_density(data = filter(recap_raw, ReleaseYear == RL), aes(x = Length,col = "Recaptures",fill = "Recaptures"), lwd = 1.4, lty =1, bw =.5, alpha = .1)+
      theme(strip.background = element_rect(fill = "white", color ="white"),
            text = element_text(size = 26),
            legend.position=c(.26,.87),
            legend.direction = "vertical",
            legend.title=element_blank(),
            legend.box.background = element_rect(fill = "transparent", color = "transparent"),
            strip.placement = "outside",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            plot.margin = unit(c(2,2,2,2), "mm"),
            panel.background = element_rect(fill = "transparent"), # bg of the panel
            plot.background = element_rect(fill = "transparent", color = NA))
  ddd <- rbind.data.frame(
    filter(d, ReleaseYear == RL) %>% dplyr::select(ReleaseYear, Length) %>% mutate(series = "Releases"),
    filter(recap_raw, ReleaseYear == RL)%>% dplyr::select(ReleaseYear, Length) %>% mutate(series = "Recaptures"))
  plength2[[i]] <- ggplot(data=ddd, aes(x = Length)) +
   stat_bin(aes(y=..density.., fill = series), binwidth  =1, position = "dodge") +
    #geom_density(aes(col = "Releases",fill = "Releases"), bw = .5, lwd = 1.4, alpha = .1) +
    #      geom_density(aes(x=Length,
    #                      y=..density..))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,55,5),
                       minor_breaks= 0:55,
                       limits = c(21.5,45.5)) +
    geom_vline(xintercept = 35, lty =2, col = "black")+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,.3,0.05), labels = paste0(seq(0,30,5), "%")) +
    theme_bw() +
    scale_color_manual(name = "", values = c("red", "orange"))+
    scale_fill_manual(name = "", values = c("red", "orange"))+
    #geom_density(data = filter(recap_raw, ReleaseYear == RL), aes(x = Length,col = "Recaptures",fill = "Recaptures"), lwd = 1.4, lty =1, bw =.5, alpha = .1)+
    theme(strip.background = element_rect(fill = "white", color ="white"),
          text = element_text(size = 26),
          legend.position=c(.26,.87),
          legend.direction = "vertical",
          legend.title=element_blank(),
          legend.box.background = element_rect(fill = "transparent", color = "transparent"),
          strip.placement = "outside",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(2,2,2,2), "mm"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA))
}
p1[[1]]
plength2[[i]]

# scale_size(name = "#Released")
#ggsave(p1,"figures/1_releases.png", device = "png", width = 12, height = 9)

# ------ RECAPTURES ----------------
#------------------------------------
recap_raw <- d %>% filter(!is.na(RecaptureDate)) %>% left_join(Catch_data, by = "CatchID") %>%
  mutate(RecaptureYear = ifelse(lubridate::month(RecaptureDate)>2, lubridate::year(RecaptureDate), lubridate::year(RecaptureDate)-1),
         YearsOut = RecaptureYear - ReleaseYear) %>% filter(YearsOut<3, between(ReleaseYear,2014,2020),
                                                            !is.na(ICES_Rectangle.y)) #%>%
  #filter(RecaptureYear <2021)

recap_raw <- mutate(recap_raw, rmonth = lubridate::month(RecaptureDate),
                    rday = ifelse(rmonth<3,365+lubridate::yday(RecaptureDate),lubridate::yday(RecaptureDate)),
                    lon = mapplots::ices.rect(ICES_Rectangle.y)$lon,
                    lat = mapplots::ices.rect(ICES_Rectangle.y)$lat
) %>% filter(!(rmonth %in% 3:6))

  #mutate(YearsOut = ifelse(YearsOut<3, YearsOut, "3+"))
# Calculate range of longitudes and latitudes
recap_raw %>% dplyr::select(ReleaseYear,RecaptureDate,RecaptureYear,YearsOut, Length)
recap <- recap_raw %>%
  dplyr::select(ICES_Rectangle.y, ReleaseYear, YearsOut) %>%
  count(ICES_Rectangle.y, ReleaseYear, YearsOut)
recap$lon <- mapplots::ices.rect(recap$ICES_Rectangle.y)$lon
recap$lat <- mapplots::ices.rect(recap$ICES_Rectangle.y)$lat
xlim = range(recap$lon, na.rm = T)
ylim = range(recap$lat, na.rm = T)

recap <- filter(recap, !is.na(ICES_Rectangle.y))



# Get map:
m <- ggplot2::map_data("world",
                       xlim = xlim + c(-.5,.5),
                       ylim = ylim +c(-.5,.5) )

dat_recap_count <- recap %>%
   group_by(ReleaseYear, YearsOut) %>% summarize(N = sum(n))


# Create plot:
dat_recap<- left_join(recap, dat_recap_count, by = c("ReleaseYear", "YearsOut"))
dat_recap <- mutate(dat_recap, YearN = paste0(ReleaseYear, "+", YearsOut, " years out, N = ", N))

dat_recap <- dat_recap %>% filter(!is.na(ICES_Rectangle.y)) %>% mutate(
  lon = mapplots::ices.rect(ICES_Rectangle.y)$lon,
  lat = mapplots::ices.rect(ICES_Rectangle.y)$lat
)
dat_recap <- dat_recap %>% mutate(ncat = factor(cut(n,
                                             breaks = c(-Inf, 1,5,10, Inf),
                                             labels = c("1","2-5","5-10", ">10"))))
p2 <- list()
timing.fig <- list()
for(i in 1:length(unique(dat_recap$ReleaseYear))){
  RY <- (2014:2020)[i]
  ddf <- filter(recap_raw, ReleaseYear == RY) %>%
    dplyr::select(rday, lon,lat) %>% group_by(lon,lat) %>% summarize(avgday = floor(mean(rday))) %>%
    mutate(rdate = as.Date(avgday, origin = "2016-01-01"),
           rmonth = lubridate::month(rdate),
           rmonth = factor(ifelse(rmonth %in% 1:2, 12+rmonth, rmonth), levels = 7:14))

  timing.fig[[i]] <- ggplot() +
    # geom_hline(yintercept = 62, col = grey(.8), lty = 2)+
    # geom_segment(data = data.frame(x = -10,xend=-10,y=62,yend=Inf),
    #              aes(x=x,xend=xend,y=y,yend=yend), col = grey(.8), lty = 2)+
    geom_polygon(data = filter(countymap, !hole), aes(long, lat, group = group),
                 fill = NA, col = "lightblue",
                 lty = 1, lwd = .4)+
    geom_polygon(data = m, aes(long, lat, group = group), fill = "grey")+
    geom_tile(data = ddf, aes(lon, lat, fill = rmonth))+
    coord_quickmap(xlim = xlim + c(0,0), ylim = ylim +c(-1,1)) +
  #  facet_grid(ReleaseYear~YearsOut) +
    theme_bw()+
  scale_fill_viridis_d("Month", direction = 1,drop=FALSE)+
    #scale_fill_manual(values = c("yellow", "orange", "red","brown3"),# colorRamps::matlab.like2(5),
    #                  name = "#Recaptures")+
    theme(strip.background = element_rect(fill = "white", color = "white"),
          #strip.text = element_text(size = 14),
          legend.position = "top",
          text = element_text(size = 26),
          panel.grid = element_blank(),
          strip.text.y = element_blank(), axis.title = element_blank(),
          plot.margin = unit(c(2,2,2,2), "mm"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA))+
    scale_x_continuous("Longitude", breaks = seq(-30,10,10),
                       labels = c("30°W","20°W", "10°W", "0°", "10°E"))+
    scale_y_continuous("Latitude", breaks = seq(50,70,5),
                       labels = paste0(seq(50,70,5), "°"))+
  #  geom_text(data= df.count, aes(label = Nlab, x = lon, y = lat), hjust = 1.05, vjust = 1.2, size =6)

  for(j in 0:2){

  k = j*length(unique(dat_recap$ReleaseYear))+i
  if(RY+j<=2020)
  {
    df <- filter(dat_recap, ReleaseYear == RY, YearsOut == j)
    df.count <- summarize(df, N = sum(n)) %>% mutate(lon = -Inf, lat = Inf, Nlab = paste0("N=",N))
  p2[[k]] <-    ggplot() +
  #     geom_hline(yintercept = 62, col = grey(.8), lty = 2)+
  # geom_segment(data = data.frame(x = -10,xend=-10,y=62,yend=Inf),
  #              aes(x=x,xend=xend,y=y,yend=yend), col = grey(.8), lty = 2)+
    geom_polygon(data = filter(countymap, !hole), aes(long, lat, group = group),
                 fill = NA, col = "lightblue",
                 lty = 1, lwd = .4)+
  geom_polygon(data = m, aes(long, lat, group = group), fill = "grey")+
  geom_tile(data = df, aes(lon, lat, fill = ncat))+
  coord_quickmap(xlim = xlim + c(0,0), ylim = ylim +c(-1,1)) +
  facet_grid(ReleaseYear~YearsOut) +
  theme_bw()+
 # guides(fill = guide_colorbar(barheight = 35))+
  # scale_fill_gradientn(colours = heat.colors(50)[45:5],
  #                      name = "#Recaptures",
  #                      limits = range(dat_recap$n)) +
  # scale_fill_manual(values = colorRamps::matlab.like2(5),
  #                   name = "#Recaptures")+
  scale_fill_manual(values = c("yellow", "orange", "red","brown3"),# colorRamps::matlab.like2(5),
                    name = "#Recaptures")+
  #scale_fill_viridis_d()+
  theme(strip.background = element_rect(fill = "white", color = "white"),
        #strip.text = element_text(size = 14),
        legend.position = "top",
        text = element_text(size = 26),
        panel.grid = element_blank(),
        strip.text.y = element_blank(), axis.title = element_blank(),
        plot.margin = unit(c(2,2,2,2), "mm"),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA))+
  scale_x_continuous("Longitude", breaks = seq(-30,10,10),
                     labels = c("30°W","20°W", "10°W", "0°", "10°E"))+
  scale_y_continuous("Latitude", breaks = seq(50,70,5),
                      labels = paste0(seq(50,70,5), "°"))+
      geom_text(data= df.count, aes(label = Nlab, x = lon, y = lat),  hjust = -.8, vjust = 1.2, size = 9)
  if(k>1) p2[[k]] <- p2[[k]]+theme(legend.position = "none")
  }else{p2[[k]] <- NULL}
  }
}
# ------------------------
plength <- plength2
library(cowplot)
legendz<- get_legend(p1[[1]]+theme(legend.text = element_text(size = 20), legend.title = element_text(size=25), legend.direction = "horizontal"))
legendz2 <- get_legend(p2[[1]]+theme(legend.text = element_text(size= 20), legend.title = element_text(size=25), legend.direction = "horizontal"))
legendz3 <- get_legend(timing.fig[[3]]+theme(legend.text = element_text(size= 20),
                                             legend.title = element_text(size=25),
                                             legend.direction = "horizontal"))
q1 <- plot_grid(
  NULL,
          p1[[1]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.x = element_text(size = 26),strip.text.y = element_text(size = 26),
                          axis.text.x = element_blank())+facet_grid("2014"~"Releases",switch = "y"),
          p1[[2]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.y = element_text(size = 26), axis.text.x = element_blank()),
          p1[[3]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.y = element_text(size = 26), axis.text.x = element_blank()),
          p1[[4]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.y = element_text(size = 26), axis.text.x = element_blank()),
          p1[[5]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.y = element_text(size = 26), axis.text.x = element_blank()),
          p1[[6]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.y = element_text(size = 26), axis.text.x = element_blank()),
          p1[[7]] + theme(plot.margin = unit(c(.1,.1,.1,.1), "mm"),legend.position = "none", strip.text.y = element_text(size = 26)),
  legendz,
          # ----
          plength[[1]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),axis.text.x = element_blank(),legend.text = element_text(size = 20), strip.text =   element_text(size = 26))+facet_grid(~"Length"),
          plength[[2]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),axis.text.x = element_blank(),legend.position = "none", strip.text.y = element_text(size = 26)),
          plength[[3]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),axis.text.x = element_blank(),legend.position = "none", strip.text.y = element_text(size = 26)),
          plength[[4]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),axis.text.x = element_blank(),legend.position = "none", strip.text.y = element_text(size = 26)),
          plength[[5]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),axis.text.x = element_blank(),legend.position = "none", strip.text.y = element_text(size = 26)),
          plength[[6]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),axis.text.x = element_blank(),legend.position = "none", strip.text.y = element_text(size = 26)),
          plength[[7]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.y = element_text(size = 26)),
 NULL,
          # --
          p2[[1]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"), legend.position = "none", strip.text.x = element_text(size = 26), axis.text.x = element_blank()) +
   facet_grid(~"Migration cycle 1"),
          p2[[2]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(), axis.text.x = element_blank()),
          p2[[3]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(), axis.text.x = element_blank()),
          p2[[4]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(), axis.text.x = element_blank()),
          p2[[5]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(), axis.text.x = element_blank()),
          p2[[6]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(), axis.text.x = element_blank()),
          p2[[7]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank()),
 legendz2,
          # ----
          p2[[8]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.x = element_text(size = 26),axis.text = element_blank())+
   facet_grid(~"Migration cycle 2"),
          p2[[9]]  + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(),axis.text = element_blank()),
          p2[[10]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(),axis.text = element_blank()),
          p2[[11]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(),axis.text = element_blank()),
          p2[[12]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(),axis.text = element_blank()),
          p2[[13]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text = element_blank(),axis.text = element_blank()),
   NULL,
 NULL,
          # ---------
          p2[[14]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),
                           legend.position = "none", strip.text.y = element_text(size = 26, angle = -90),
                           strip.text.x = element_text(size = 26),axis.text = element_blank())+
   facet_grid(~"Migration cycle 3"),
          p2[[15]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.x = element_blank(),axis.text = element_blank()),
          p2[[16]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.x = element_blank(),axis.text = element_blank()),
          p2[[17]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.x = element_blank(),axis.text = element_blank()),
          p2[[18]] + theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.x = element_blank(),axis.text = element_blank()),
#          p2[[19]] + theme(legend.position = "none"),
        #  legendz2,
        #  legendz3,
NULL,NULL,legendz3,
timing.fig[[1]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", strip.text.x = element_text(size = 26),
                      axis.text = element_blank()) +
  facet_grid(~"Timing"),
timing.fig[[2]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", axis.text = element_blank()),
timing.fig[[3]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", axis.text = element_blank()),
timing.fig[[4]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", axis.text = element_blank()),
timing.fig[[5]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", axis.text = element_blank()),
timing.fig[[6]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", axis.text = element_blank()),
timing.fig[[7]]+theme(plot.margin = unit(c(.5,0,.5,0), "mm"),legend.position = "none", axis.text.y = element_blank()),

#plot_grid(legendz, legendz2,ncol = 2),
          ncol = 6, nrow = 8, byrow =F,
          rel_heights = c(.25,1.14,rep(1,5),1.08),
            rel_widths = c(.62,.9,.93,rep(.85,3)))
ggsave(q1, file = "maps_of_releases_and_recaptures_2.pdf", width =7900, height = 8000, unit = "px")


# # ---------------------------
#
#  # scale_size(name = "#Released")
# #ggsave(p2, "figures/4a_recaptures_2011-15.png", device = "png",  width = 9, height = 9)
#
#
# # ------------------------------------
# ggplot() +
#   geom_hline(yintercept = 62, col = grey(.8), lty = 2)+
#   geom_segment(data = data.frame(x = -10,xend=-10,y=62,yend=Inf),
#                aes(x=x,xend=xend,y=y,yend=yend), col = grey(.8), lty = 2)+
#   geom_polygon(data = m, aes(long, lat, group = group), fill = "grey")+
#   geom_tile(data = dat_recap, aes(lon, lat, fill = ReleaseMonth))+
#   coord_quickmap(xlim = xlim + c(0,0), ylim = ylim +c(-1,1)) +
#   facet_grid(~ReleaseYear) +
#   theme_bw()+
#   # guides(fill = guide_colorbar(barheight = 35))+
#   # scale_fill_gradientn(colours = heat.colors(50)[45:5],
#   #                      name = "#Recaptures",
#   #                      limits = range(dat_recap$n)) +
#   # scale_fill_manual(values = colorRamps::matlab.like2(5),
#   #                   name = "#Recaptures")+
#   scale_fill_manual(values = c("yellow", "orange", "red","brown3"),# colorRamps::matlab.like2(5),
#                     name = "#Recaptures")+
#   #scale_fill_viridis_d()+
#   theme(strip.background = element_rect(fill = "white", color = "white"),
#         #strip.text = element_text(size = 14),
#         legend.position = "top",
#         text = element_text(size = 18),
#         panel.grid = element_blank(),
#         strip.text.y = element_blank(), axis.title = element_blank())+
#   scale_x_continuous("Longitude", breaks = seq(-30,10,10),
#                      labels = c("30°W","20°W", "10°W", "0°", "10°E"))+
#   scale_y_continuous("Latitude", breaks = seq(50,70,5),
#                      labels = paste0(seq(50,70,5), "°N"))
#
#
# # ----------------------------
#
#
#
#
#
# library(ggpubr)
# ggarrange(p1+ theme(legend.position = "left"),
#           p2, ncol = 2)
#
# # -----------------------------
# recap <- d %>% filter(!is.na(RecaptureDate)) %>% left_join(Catch_data, by = "CatchID") %>%
#   mutate(RecaptureYear = lubridate::year(RecaptureDate),
#          YearsOut = RecaptureYear - ReleaseYear) %>%
#   filter(RecaptureYear <2021,
#          !is.na(ICES_Rectangle.y))
#
# d.count <- recap %>%
#   group_by(RecaptureYear, ICES_Rectangle.y) %>% count()
# d.count$lon <- mapplots::ices.rect(d.count$ICES_Rectangle.y)$lon
# d.count$lat <- mapplots::ices.rect(d.count$ICES_Rectangle.y)$lat
# #d.count <- filter(d.count, lat<75)
#
# xlim = range(d.count$lon, na.rm = T)
# ylim = range(d.count$lat, na.rm = T)
# m <- ggplot2::map_data("world",
#                        xlim = xlim + c(-.5,.5),
#                        ylim = ylim +c(-.5,.5) )
#
# recap_count <- d.count %>%
#   group_by(RecaptureYear) %>% summarize(N = sum(n))
#
#
# # Create plot:
# d.count<- left_join(d.count, recap_count, by = "RecaptureYear")
# d.count <- mutate(d.count, YearN = paste0(RecaptureYear, ", ", N, " recaptures"))
# d.count <- mutate(d.count,
#                   catch.cat = cut(n,
#                                   breaks = c(-Inf, 1, 5, 10, Inf),
#                                   labels=c("1","2-5","6-10", ">10")))
#
# ggplot() +
#   geom_polygon(data = m, aes(long, lat, group = group), fill = "grey") +
#   geom_tile(data = d.count, aes(lon, lat, fill = catch.cat))+#,# col = CatchWeight)) +
#   coord_quickmap(xlim = xlim + c(-.5,.5), ylim = ylim +c(-.5,.5)) +
#   facet_wrap( ~YearN, ncol = 3) +
#   theme_bw()+
#   # guides(fill = guide_colorbar(barheight = 55))+
#   scale_fill_manual(values = c("yellow", "orange", "red", "brown3"), #heat.colors(4)[4:1],
#                     name = "Recaptures") +
#   theme(strip.background = element_rect(fill = "white"),
#         text = element_text(size = 14),
#         legend.position = "top")+
#   scale_x_continuous("Longitude", breaks = seq(-30,10,10),
#                      labels = c("30?W","20?W","10?W", "0?", "10?E"))+
#   scale_y_continuous("Latitude", breaks = seq(50,90,5),
#                      labels = paste0(seq(50,90,5), "?N"))
# ggsave("figures/3_Recaptures_by_year.png", device = "png", width =9, height = 9, dpi = "retina")
#
#
# # -----------------------------------------
# # Catch data
# recap <- d %>% filter(!is.na(RecaptureDate)) %>% left_join(Catch_data, by = "CatchID") %>%
#   mutate(RecaptureYear = lubridate::year(RecaptureDate),
#          YearsOut = RecaptureYear - ReleaseYear) %>%
#   filter(between(RecaptureYear,2014,2020),
#          !is.na(ICES_Rectangle.y))
# recap$lon <- mapplots::ices.rect(recap$ICES_Rectangle.y)$lon
# recap$lat <- mapplots::ices.rect(recap$ICES_Rectangle.y)$lat
#
# ggplot() +
#   geom_polygon(data = m, aes(long, lat, group = group), fill = "grey") +
#   geom_jitter(data = recap, aes(lon, lat, col = factor(lubridate::month(RecaptureDate))))+#,# col = CatchWeight)) +
#   coord_quickmap(xlim = xlim + c(-.5,.5), ylim = ylim +c(-.5,.5)) +
#   facet_wrap( ~RecaptureYear,nrow = 2) +
#   theme_bw()+
#   # guides(fill = guide_colorbar(barheight = 55))+
#   #scale_fill_manual(values = colorspace::heat_hcl(12))+
#   scale_color_viridis_d(direction = -1, name = "Month")+
#   theme(strip.background = element_rect(fill = "white"),
#         text = element_text(size = 14),
#         legend.position = "top",
#         panel.grid = element_blank())+
#   scale_x_continuous("Longitude", breaks = seq(-30,10,10),
#                      labels = c("30°W","20°W","10°W", "0°", "10°E"))+
#   scale_y_continuous("Latitude", breaks = seq(50,90,5),
#                      labels = paste0(seq(50,90,5), "°N"))
# #---
#
#
#
#
#
#
# d.count <- d %>% mutate(CatchYear = lubridate::year(CatchDate)) %>%
#   group_by(CatchYear, ICES_Rectangle) %>% summarize(n = sum(CatchWeight))
# d.count$lon <- mapplots::ices.rect(d.count$ICES_Rectangle)$lon
# d.count$lat <- mapplots::ices.rect(d.count$ICES_Rectangle)$lat
# #d.count <- filter(d.count, lat<75)
#
# xlim = range(d.count$lon, na.rm = T)
# ylim = range(d.count$lat, na.rm = T)
# m <- ggplot2::map_data("world",
#                        xlim = xlim + c(-.5,.5),
#                        ylim = ylim +c(-.5,.5) )
#
# catch_count <- d.count %>%
#   group_by(CatchYear) %>% summarize(N = sum(n))
#
#
# # Create plot:
# d.count<- left_join(d.count, catch_count, by = "CatchYear")
# d.count <- mutate(d.count, YearN = paste0(CatchYear, ", ", round(N/1e3,0), " t"))
# d.count <- mutate(d.count,
#                   catch.cat = cut(n/1e3,
#                                   breaks = c(-Inf, 1000, 3000, 5000, Inf),
#                                   labels=c("<1000","1000-3000","3000-5000", ">5000"))) %>%
#   filter(CatchYear <=2020)
#
# ggplot() +
#   geom_polygon(data = m, aes(long, lat, group = group), fill = "grey") +
#   geom_tile(data = d.count, aes(lon, lat, fill = catch.cat))+#,# col = CatchWeight)) +
#   coord_quickmap(xlim = xlim + c(-.5,.5), ylim = ylim +c(-.5,.5)) +
#   facet_wrap( ~YearN, ncol = 3) +
#   theme_bw()+
#  # guides(fill = guide_colorbar(barheight = 55))+
#   scale_fill_manual(values = c("yellow", "orange", "red", "brown3"), #heat.colors(4)[4:1],
#                        name = "Catch (t)") +
#   theme(strip.background = element_rect(fill = "white"),
#         text = element_text(size = 14),
#         legend.position = "top",
#         panel.grid = element_blank())+
#   scale_x_continuous("Longitude", breaks = seq(-30,10,10),
#                      labels = c("30°W","20°W","10°W", "0°", "10°E"))+
#   scale_y_continuous("Latitude", breaks = seq(50,90,5),
#                      labels = paste0(seq(50,90,5), "°N"))
# ggsave("2_CatchData.png", device = "png", width = 10, height = 9, dpi = "retina")
#

# ------------------------------
recap <- recap_raw %>%
  dplyr::select(ICES_Rectangle.y, ReleaseYear, YearsOut, Length, Latitude) #%>%
  #count(ICES_Rectangle.y, ReleaseYear, YearsOut, Length)
recap$lon <- mapplots::ices.rect(recap$ICES_Rectangle.y)$lon
recap$lat <- mapplots::ices.rect(recap$ICES_Rectangle.y)$lat
xlim = range(recap$lon, na.rm = T)
ylim = range(recap$lat, na.rm = T)

recap <- filter(recap, !is.na(ICES_Rectangle.y))



# Get map:
m <- ggplot2::map_data("world",
                       xlim = xlim + c(-.5,.5),
                       ylim = ylim +c(-.5,.5) )
as.data.frame(ICESecoregion)

table(as.character(countymap$group))

ggplot(filter(countymap, group != "Greater North Sea.3242")) +
  geom_polygon(aes(long, lat, group = group), fill = NA, col = "lightblue",
                                     lty = 1, lwd = .3)+
  coord_quickmap(xlim = xlim + c(0,0), ylim = ylim +c(-1,1))


dat_recap <- recap %>% group_by(lon,lat,ICES_Rectangle.y, ReleaseYear, YearsOut) %>% summarize(
  mlat = mean(Latitude),
  mlength = mean(Length),
  n = n())
# Create plot:
#dat_recap <- mutate(recap, YearN = paste0(ReleaseYear, "+", YearsOut, " years out, N = ", N))

dat_recap <- dat_recap %>% filter(!is.na(ICES_Rectangle.y)) %>% mutate(
  lon = mapplots::ices.rect(ICES_Rectangle.y)$lon,
  lat = mapplots::ices.rect(ICES_Rectangle.y)$lat
)
dat_recap <- dat_recap %>% mutate(ncat = factor(cut(n,
                                                    breaks = c(-Inf, 1,5,10, Inf),
                                                    labels = c("1","2-5","5-10", ">10"))))
mlengthfig <- list()

  dat_recap <- dat_recap %>% mutate(cycle = paste0("Migration cycle ",YearsOut+1))

  p1 <- ggplot() +
    #geom_hline(yintercept = 62, col = grey(.8), lty = 2)+
    #geom_segment(data = data.frame(x = -10,xend=-10,y=62,yend=Inf),
    #             aes(x=x,xend=xend,y=y,yend=yend), col = grey(.8), lty = 2)+
    geom_polygon(data = filter(countymap, !hole), aes(long, lat, group = group), fill = NA, col = "lightblue",
                 lty = 1, lwd = .1)+
    geom_polygon(data = m, aes(long, lat, group = group), fill = "grey")+
    # geom_sf(data=ICESecoregion, aes(geometry = geometry, color=Ecoregion), fill = NA)+
    geom_tile(data = dat_recap, aes(lon, lat, fill = mlength))+
    coord_quickmap(xlim = xlim + c(0,0), ylim = ylim +c(-1,1)) +
    guides(fill = guide_colorbar(barheight = 25))+
    facet_grid(ReleaseYear~cycle)+
    #  facet_grid(ReleaseYear~YearsOut) +
    theme_bw()+

    scale_fill_viridis_c("Mean\nLength", direction = -1, breaks = seq(20,50,2), minor_breaks = 25:45)+
    #scale_fill_manual(values = c("yellow", "orange", "red","brown3"),# colorRamps::matlab.like2(5),
    #                  name = "#Recaptures")+
    theme(strip.background = element_rect(fill = "transparent", color = "transparent"),
          #strip.text = element_text(size = 14),
          legend.position = "right",
          #text = element_text(size = 26),
          panel.grid = element_blank(),
          axis.text = element_text(size = 8),
         # strip.text.y = element_blank(), axis.title = element_blank(),
         # plot.margin = unit(c(2,2,2,2), "mm"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA))+
    scale_x_continuous("Longitude", breaks = seq(-30,10,10),
                       labels = c("30°W","20°W", "10°W", "0°", "10°E"))+
    scale_y_continuous("Latitude", breaks = seq(50,70,5),
                       labels = paste0(seq(50,70,5), "°"))
  ggsave(p1,file="MEANLENGTH.pdf", width = 6, height = 9)
  ggplot() +
    #geom_hline(yintercept = 62, col = grey(.8), lty = 2)+
    #geom_segment(data = data.frame(x = -10,xend=-10,y=62,yend=Inf),
    #             aes(x=x,xend=xend,y=y,yend=yend), col = grey(.8), lty = 2)+
    geom_polygon(data = filter(countymap, !hole), aes(long, lat, group = group),
                 fill = NA, col = "lightblue",
                 lty = 1, lwd = .1)+
    geom_polygon(data = m, aes(long, lat, group = group), fill = "grey")+
    geom_tile(data = dat_recap, aes(lon, lat, fill = mlat))+
    coord_quickmap(xlim = xlim + c(0,0), ylim = ylim +c(-1,1)) +
    guides(fill = guide_colorbar(barheight = 25))+
    facet_grid(ReleaseYear~cycle)+
    #  facet_grid(ReleaseYear~YearsOut) +
    theme_bw()+
    scale_fill_viridis_c("Mean\nRelease\nLatitude", direction = -1)+#, breaks = seq(20,50,2), minor_breaks = 25:45)+
    #scale_fill_manual(values = c("yellow", "orange", "red","brown3"),# colorRamps::matlab.like2(5),
    #                  name = "#Recaptures")+
    theme(strip.background = element_rect(fill = "white", color = "white"),
          #strip.text = element_text(size = 14),
          legend.position = "right",
          #text = element_text(size = 26),
          panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          # strip.text.y = element_blank(), axis.title = element_blank(),
          # plot.margin = unit(c(2,2,2,2), "mm"),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA))+
    scale_x_continuous("Longitude", breaks = seq(-30,10,10),
                       labels = c("30°W","20°W", "10°W", "0°", "10°E"))+
    scale_y_continuous("Latitude", breaks = seq(50,70,5),
                       labels = paste0(seq(50,70,5), "°"))
  ggsave("MEANLATITUDE.pdf", width = 6, height = 9)
    #  geom_text(data= df.count, aes(label = Nlab, x = lon, y = lat), hjust = 1.05, vjust = 1.2, size =6)
