map_area <- get_map(location= c(left   = min(Data_mackerel_final$cLon, Data_mackerel_final$lo, na.rm=TRUE)-1,
                                bottom = min(Data_mackerel_final$cLat, Data_mackerel_final$la, na.rm=TRUE)-1,
                                right  = max(Data_mackerel_final$cLon, Data_mackerel_final$lo, na.rm=TRUE)+1,
                                top    = max(Data_mackerel_final$cLat, Data_mackerel_final$la, na.rm=TRUE)+1))

plot_timerelease_simple <- function(Year)
{
  lag = 1
  dat_release <- Data_mackerel_all %>% subset(Release_year==Year &
                                              Tag_area %in% c("West_Ireland", "North_Ireland"))
dat_release$julian <-  as.numeric(julian(dat_release$relesedate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
dat_release$Tag_area <- as.character(dat_release$Tag_area)
dat_release$Tag_area <- as.factor(dat_release$Tag_area)
dat_release$Tag_area <- factor(dat_release$Tag_area, labels=c("From N. Ireland", "From W. Ireland"))
dat_recap <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                              Tag_area %in% c("West_Ireland", "North_Ireland"))
dat_recap$julian <-  as.numeric(julian(dat_recap$recapturedate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
dat_recap$Tag_area <- as.character(dat_recap$Tag_area)
dat_recap$Tag_area <- as.factor(dat_recap$Tag_area)

dat_recap_new <- as_tibble(dat_recap) %>% dplyr::select(cLon, cLat, julian) %>% count(round(cLon,1), round(cLat,1), round(julian,0))
colnames(dat_recap_new) <- c("cLon", "cLat", "julian_recapture", "n")

# lag 2
if(Year < 2018){
lag = 2
dat_recap2 <- Data_mackerel_final %>% subset(Release_year==Year & Catch_year==Year+lag &
                                              Tag_area %in% c("West_Ireland", "North_Ireland"))
dat_recap2$julian <-  as.numeric(julian(dat_recap2$recapturedate, as.POSIXct(paste0(Year, "-01-01"), tz = "GMT")))
dat_recap2$Tag_area <- as.character(dat_recap2$Tag_area)
dat_recap2$Tag_area <- as.factor(dat_recap2$Tag_area)

dat_recap2_new <- as_tibble(dat_recap2) %>% dplyr::select(cLon, cLat, julian) %>% count(round(cLon,1), round(cLat,1), round(julian,0))
colnames(dat_recap2_new) <- c("cLon", "cLat", "julian_recapture", "n")


dat_release_new <- as_tibble(dat_release) %>% dplyr::select(lo, la, julian) %>% count(round(lo,1), round(la,1), round(julian,0))
colnames(dat_release_new) <- c("lo", "la", "julian_release", "n")
ncount_recap2 <- dat_recap2 %>% count()
ncount_recap2$lon = c(-10)
ncount_recap2$lat = c(70)
ncount_recap2$label = paste0("n=", ncount_recap2$n)
}

ncount_release <- dat_release %>% count()
ncount_release$lon = c(-32)
ncount_release$lat = c(53)
ncount_release$label = paste0("#Released: ", ncount_release$n, "\\#Recaptured: ",ncount_recap$n )
ncount_recap <- dat_recap %>% count()
ncount_recap$lon = c(-10)
ncount_recap$lat = c(70)
ncount_recap$label = paste0("n=", ncount_recap$n)

if(FALSE){
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
}
#dat_recap_new<-dat_recap_new %>% rename(n =N)

if(Year < 2018){
plot_data <- data.frame(
  lo = c(dat_release_new$lo, dat_recap_new$cLon, dat_recap2_new$cLon),
  la = c(dat_release_new$la, dat_recap_new$cLat,dat_recap2_new$cLat),
  julian = c(dat_release_new$julian_release, dat_recap_new$julian_recapture, dat_recap2_new$julian_recapture),
  n = c(dat_release_new$n, dat_recap_new$n, dat_recap2_new$n),
  type = c(rep("Released", nrow(dat_release_new)), rep("Recaptured year 1 after release", nrow(dat_recap_new)),
           rep("Recaptured year 2 after release", nrow(dat_recap2_new)))
)


g1 <- ggmap(map_area) +
  geom_point(data = plot_data[plot_data$type == "Released",], aes(x=lo, y=la, size=n),col = "blue", shape=19)+
  geom_point(data = plot_data[plot_data$type == "Recaptured year 1 after release",],aes(x=lo, y=la, col = type), size = 2)+
  geom_point(data = plot_data[plot_data$type == "Recaptured year 2 after release",],aes(x=lo, y=la, col = type), size = 1)+

  #,col=julian_release)) +
  #geom_point(data=plot_data, aes(x=lo, y=la, size=n, col = type), shape=19)
  #geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat)) +
 # scale_fill_viridis_c() +
  scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6), name = "#Released") +
  scale_color_manual(name = "",values = c("red", "green", "blue"))+
  guides(size =  guide_legend(override.aes = list(col = "blue", fill=NA)))+
  theme_bw()+theme(legend.position = "bottom")+
  #new_scale_color() +
 # geom_point(data=dat_recap_new, aes(x=cLon, y=cLat),#, fill=julian_recapture),
 #            size = 2, col = "red") +
  #scale_fill_viridis_c() +
  #scale_size_continuous(range = c(0.3, 0.6)) +
  #geom_text(aes(label = paste(str_pad("#Released:",width = 19, side = "right", " "), ncount_release$n,"\n",
                            #"#Recaptured ",Year+1, ":",str_pad(ncount_recap$n,side = "left",width = 4, " "),"\n",
                            #"#Recaptured ",Year+2, ":",str_pad(ncount_recap2$n,side = "left",width = 4, " "),"\n",
                            #sep = ""),
               # x = Inf, y = Inf), hjust =1, vjust = 1.1, col = "black")+
 # geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"),
  #          aes(x=x, y=y, label=label), col="black", fontface="bold") +
 # geom_text(data=data.frame(x=-26,y=55,label="Nsamp release:"),
  #          aes(x=x, y=y, label=label), col="black", fontface="bold") +
  ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5))
}
if(Year == 2018){
  plot_data <- data.frame(
    lo = c(dat_release_new$lo, dat_recap_new$cLon),
    la = c(dat_release_new$la, dat_recap_new$cLat),
    julian = c(dat_release_new$julian_release, dat_recap_new$julian_recapture),
    n = c(dat_release_new$n, dat_recap_new$n),
    type = c(rep("Released", nrow(dat_release_new)), rep("Recaptured year 1 after release", nrow(dat_recap_new)))
  )


  g1 <- ggmap(map_area) +
    geom_point(data=plot_data, aes(x=lo, y=la, size=n, col = type), shape=19)+
    geom_point(data = plot_data[plot_data$type == "Recaptured year 1 after release",],aes(x=lo, y=la, col = type), size = 2)+
    #geom_point(data = plot_data[plot_data$type == "Recaptured lag 2",],aes(x=lo, y=la, col = type), size = 1)+
    #,col=julian_release)) +
    #geom_point(data=plot_data, aes(x=lo, y=la, size=n, col = type), shape=19)
    #geom_polygon(data = hull, alpha = 0.5, aes(x=cLon, y=cLat)) +
    # scale_fill_viridis_c() +
    scale_size_continuous(breaks=c(1,10,100,1000),range = c(0.3, 6), name = "#Released") +
    scale_color_manual(name = "",values = c("red", "blue"))+
    guides(size =  guide_legend(override.aes = list(col = "blue", fill=NA)))+
    theme_bw()+theme(legend.position = "bottom")+
    #new_scale_color() +
    # geom_point(data=dat_recap_new, aes(x=cLon, y=cLat),#, fill=julian_recapture),
    #            size = 2, col = "red") +
    #scale_fill_viridis_c() +
    #scale_size_continuous(range = c(0.3, 0.6)) +
   # geom_text(aes(label = paste(str_pad("#Released:",width = 19, side = "right", " "), ncount_release$n,"\n",
    #                            "#Recaptured ",Year+1, ":",str_pad(ncount_recap$n,side = "left",width = 4, " "),"\n",
    #                            #"#Recaptured ",Year+1, ":",str_pad(ncount_recap2$n,side = "left",width = 3, " "),"\n",
    #                            sep = ""),
    #              x = Inf, y = Inf), hjust =1, vjust = 1.1, col = "black")+
    # geom_text(data=data.frame(x=-8,y=71,label="Nsamp recapture:"),
    #          aes(x=x, y=y, label=label), col="black", fontface="bold") +
    # geom_text(data=data.frame(x=-26,y=55,label="Nsamp release:"),
    #          aes(x=x, y=y, label=label), col="black", fontface="bold") +
    ggtitle(Year) + theme(plot.title = element_text(hjust = 0.5))

}

g1

RES <- list()
RES$plot <- g1
RES$area <- area
return(RES)
}


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

Years = c(2011:2018)
SAVE <- list()
  # areas <- c()
for (yr in 1:length(Years)){
  SAVE[[yr]] <- plot_timerelease_simple(Years[[yr]])
  # areas <- rbind(areas, SAVE[[yr]][[2]])
}
my.legend <- g_legend(SAVE[[1]]$plot+theme(legend.text = element_text(size = 12),
                                           legend.title= element_text(size = 12)))
  # areas <- c()
for (yr in 1:length(SAVE)){
  SAVE[[yr]]$plot <- SAVE[[yr]]$plot+theme(legend.position = "none",
                                           plot.margin = unit(c(0,0,0,2), "mm"))+xlab("")+ylab("")

}
ggg_timerelease_size1 <- grid.arrange(my.legend,arrangeGrob(grobs=map(SAVE, pluck, "plot"),
                                                             byrow =FALSE,ncol=2, nrow = 4,
                                                             layout_matrix = matrix(1:8, ncol = 2, nrow= 4, byrow=TRUE),
                                                             padding = unit(0, "line")),
                                        nrow =2, heights = c(4,100))
ggsave(ggg_timerelease_size1, filename=paste0("../plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag12", "_part1.png"),
       width=25, height=40, units="cm", device = "png", dpi = 400)
ggsave(ggg_timerelease_size1, filename=paste0("../plots/mark_recap_raw_", Years[1], "_", tail(Years, 1), "_lag12", "_part1.pdf"),
       width=25, height=40, units="cm", device = "pdf", dpi = 400)


