##' Function to create the marginal effect plot
##' @param model is the name of the saved model object from which to extract the model run and data
##' @param area_sel is the "go-to" direction to choose: choice among the four ICES ecoregions
##' @param var  variable name to match the one in the data frame
##' @param varlabel  the variable label to show in the plot
##' @param ybreaks  the breaks for the y axis. needed for plotting size adjustement
##' @param ylabels  the labels for the y axis. needed for plotting size adjustement
##' @param fix_date  to you want to fix the date for the marginal effect plot? if so specify in mean_date the date to use
##' @param mean_date  the mean date to use for the marginal effect plot
##' @param center  center the variable y-y[1]
##' @param last  is this the last row of figure? controls plotting output
##' @details creates all the individual panel that aggregates results across migration cycles
##' @return a ggplot object to be aggregated with ggarrange for final plot
##' @export create_panel
##'
create_panel <- function(model=list(Dec_0lag[[1]],Dec_1lag[[1]],Dec_2lag[[1]]),
                         area_sel='toiceland', var="Latitude", varlabel="Latitude (º)",
                         ybreaks = waiver(), fix_date = TRUE, mean_date = 60, center = TRUE,
                         ylabels = waiver(), last=FALSE){
  which.col <- which(colnames(model[[1]]$data) == area_sel)
  # migration cycle 1
  ICEdat <- model[[1]]$data[which(model[[1]]$data[,which.col] == 1),]
  vals = ifelse(fix_date == TRUE, (mean_date-mean(model[[1]]$data$julian_recapture_std))/sd(model[[1]]$data$julian_recapture_std), mean(ICEdat$julian_recapture_std_scaled))
  asd <- visreg::visreg(fit=model[[1]]$m1, xvar=var, plot=FALSE, type = "conditional",
                 data=model[[1]]$data, scale="linear",
                 cond = list(Latitude = 52.5, #mean(ICEdat$Latitude),
                             Release_year = 2017,
                             Length= 35, #mean(ICEdat$Length),
                             julian_recapture_std_scaled = vals))
  if (center== TRUE) {
    asd$fit = asd$fit %>% mutate(visregLwr = visregLwr-visregFit[1],
                                 visregUpr = visregUpr-visregFit[1],
                                 visregFit = visregFit-visregFit[1])
  }
  # migration cycle 2
  which.col1 <- which(colnames(model[[2]]$data) == area_sel)
  ICEdat1 <- model[[2]]$data[which(model[[2]]$data[,which.col1] == 1),]
  vals1 = ifelse(fix_date == TRUE, (mean_date-mean(model[[2]]$data$julian_recapture_std))/sd(model[[2]]$data$julian_recapture_std), mean(ICEdat1$julian_recapture_std_scaled))
  asd1 <- visreg::visreg(fit=model[[2]]$m1, xvar=var, plot=FALSE, type = "conditional",
                 data=model[[2]]$data, scale="linear",
                 cond = list(Latitude = 52.5, #mean(ICEdat1$Latitude),
                             Release_year = 2017,
                             Length= 35, #mean(ICEdat1$Length),
                             julian_recapture_std_scaled = vals1))
  if (center== TRUE) {
    asd1$fit = asd1$fit %>% mutate(visregLwr = visregLwr-visregFit[1],
                                   visregUpr = visregUpr-visregFit[1],
                                   visregFit = visregFit-visregFit[1])
  }
  # migration cycle 3
  which.col2 <- which(colnames(model[[3]]$data) == area_sel)
  ICEdat2 <- model[[3]]$data[which(model[[3]]$data[,which.col2] == 1),]
  vals2 = ifelse(fix_date == TRUE, (mean_date-mean(model[[3]]$data$julian_recapture_std))/sd(model[[3]]$data$julian_recapture_std), mean(ICEdat2$julian_recapture_std_scaled))
  asd2 <- visreg::visreg(fit=model[[3]]$m1, xvar=var, plot=FALSE, type = "conditional",
                 data=model[[3]]$data, scale="linear",
                 cond = list(Latitude = 52.5, #mean(ICEdat2$Latitude),
                             Release_year = 2017,
                             Length= 35, #mean(ICEdat2$Length),
                             julian_recapture_std_scaled = vals2))
  if (center== TRUE) {
    asd2$fit = asd2$fit %>% mutate(visregLwr = visregLwr-visregFit[1],
                                   visregUpr = visregUpr-visregFit[1],
                                   visregFit = visregFit-visregFit[1])
  }

    if (var %in% c("Latitude", "Length")){
      my_col <- sym(var)
      out <- ggplot(data=asd$fit, aes(x=!!my_col, y=visregFit)) +
      geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill="red", alpha=0.2) +
      geom_line(size=1.5, col="red") +
      geom_line(data=asd1$fit, aes(x=!!my_col, y=visregFit), col="orange", size=1.5) +
      geom_line(data=asd2$fit, aes(x=!!my_col, y=visregFit), col="#FEE724", size=1.5) +
      theme_bw() +
      labs(y="", x=varlabel) +
      scale_y_continuous(breaks = ybreaks, labels = ylabels) +
      theme(axis.title = element_text(size=15),
            axis.text = element_text(size=13),
            axis.text.y = element_text(angle = 90, hjust=0.5),
            strip.text.x = element_text(size = 13))
    }

    if (var == "julian_recapture_std_scaled"){
      asd$fit <- asd$fit %>% mutate(xtransf = julian_recapture_std_scaled*sd(model[[1]]$data$julian_recapture_std)+mean(model[[1]]$data$julian_recapture_std))
      asd1$fit <- asd1$fit %>% mutate(xtransf = julian_recapture_std_scaled*sd(model[[2]]$data$julian_recapture_std)+mean(model[[2]]$data$julian_recapture_std))
      asd2$fit <- asd2$fit %>% mutate(xtransf = julian_recapture_std_scaled*sd(model[[3]]$data$julian_recapture_std)+mean(model[[3]]$data$julian_recapture_std))
      out <- ggplot(data=asd$fit, aes(x=xtransf, y=visregFit)) +
      geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill="red", alpha=0.2) +
      geom_line(size=1.5, col="red") +
      geom_line(data=asd1$fit, aes(x=xtransf, y=visregFit), col="orange", size=1.5) +
      geom_line(data=asd2$fit, aes(x=xtransf, y=visregFit), col="#FEE724", size=1.5) +
      theme_bw() +
      scale_y_continuous(breaks = ybreaks, labels = ylabels) +
      scale_x_continuous(name = "Recapture date", breaks=c(1,61,122,183,245), labels=c("Jun 1", "Aug 1", "Oct 1", "Dec 1", "Feb 1")) +
      labs(y="", x=varlabel) +
      theme(axis.title = element_text(size=15),
            axis.text = element_text(size=13),
            axis.text.y = element_text(angle = 90, hjust=0.5),
            strip.text.x = element_text(size = 13))
    }
    if (var == "Release_year"){
      asd$fit <- asd$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
      asd1$fit <- asd1$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
      asd2$fit <- asd2$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
      out <- ggplot(asd$fit, aes(x=Release_year_fct, y=visregFit)) + geom_point(col="red") +
        geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col="red", width=0.4, alpha=0.8) +
        geom_point(data=asd1$fit, aes(x=Release_year_fct, y=visregFit), col="orange", size=1.5) +
        geom_point(data=asd2$fit, aes(x=Release_year_fct, y=visregFit), col="#FEE724", size=1.5) +
        stat_summary(data=asd$fit, fun=sum, group=1, geom="line", size=1.5, col="red") +
        stat_summary(data=asd1$fit, fun=sum, group=1, geom="line", col="orange", size=1.5) +
        stat_summary(data=asd2$fit, fun=sum, group=1, geom="line", col="#FEE724", size=1.5) +
        theme_bw() +scale_x_discrete(drop=FALSE) +
        labs(y="", x=varlabel) +
        scale_y_continuous(breaks = ybreaks, labels = ylabels) +
        theme(axis.title = element_text(size=15),
              axis.text = element_text(size=13),
              axis.text.y = element_text(angle = 90, hjust = 0.5),
              strip.text.x = element_text(size = 13))
    }

  if (last == FALSE) {
    out <- out + theme(panel.grid.major.x = element_blank(),
                       panel.grid.major.y = element_line(linetype = "dashed"),
                       panel.grid.minor.y = element_blank(),
                       axis.title.x = element_blank(),
                       strip.background = element_rect(fill="white", color="white"))
  }
  if (last == TRUE) {
    out <- out + theme(panel.grid.major.x = element_blank(),
                       panel.grid.major.y = element_line(linetype = "dashed"),
                       panel.grid.minor.y = element_blank(),
                       strip.background = element_rect(fill="white", color="white"))
  }
  return(out)
}



pp1    <- create_panel(model=list(Dec_0lag[[1]],Dec_1lag[[1]],Dec_2lag[[1]]), ybreaks = waiver(),
                       area_sel='toiceland', var="Latitude", varlabel="Latitude (º)",
                       fix_date = TRUE, mean_date = 60)
pp2    <- create_panel(model=list(Dec_0lag[[1]],Dec_1lag[[1]],Dec_2lag[[1]]), ybreaks = waiver(),
                       area_sel='toiceland', var="Length", varlabel="Body length (cm)",
                       fix_date = TRUE, mean_date = 60)
pp3    <- create_panel(model=list(Dec_0lag[[1]],Dec_1lag[[1]],Dec_2lag[[1]]), ybreaks = waiver(),
                       area_sel='toiceland', var="julian_recapture_std_scaled", varlabel="Recapture date (std)",
                       fix_date = TRUE, mean_date = 60)
pp4    <- create_panel(model=list(Dec_0lag[[1]],Dec_1lag[[1]],Dec_2lag[[1]]), ybreaks=waiver(),
                       area_sel='toiceland', var="Release_year", varlabel="Release year",
                       fix_date = TRUE, mean_date = 60)
pp1_NO <- create_panel(model=list(Dec_0lag[[2]],Dec_1lag[[2]],Dec_2lag[[2]]), ybreaks = waiver(),
                       area_sel='to_norway', var="Latitude", varlabel="Latitude (º)",
                       fix_date = TRUE, mean_date = 90)
pp2_NO <- create_panel(model=list(Dec_0lag[[2]],Dec_1lag[[2]],Dec_2lag[[2]]), ybreaks = waiver(),
                       area_sel='to_norway', var="Length", varlabel="Body length (cm)",
                       fix_date = TRUE, mean_date = 90)
pp3_NO <- create_panel(model=list(Dec_0lag[[2]],Dec_1lag[[2]],Dec_2lag[[2]]), ybreaks = waiver(),
                       area_sel='to_norway', var="julian_recapture_std_scaled", varlabel="Recapture date (std)",
                       fix_date = TRUE, mean_date = 90)
pp4_NO <- create_panel(model=list(Dec_0lag[[2]],Dec_1lag[[2]],Dec_2lag[[2]]), ybreaks = waiver(), ylabels=waiver(),
                       area_sel='to_norway', var="Release_year", varlabel="Release year",
                       fix_date = TRUE, mean_date = 90)
pp1_NS <- create_panel(model=list(Dec_0lag[[3]],Dec_1lag[[3]],Dec_2lag[[3]]), ybreaks = waiver(),
                       area_sel='to_northsea', var="Latitude", varlabel="Latitude (º)",
                       fix_date = TRUE, mean_date = 150)
pp2_NS <- create_panel(model=list(Dec_0lag[[3]],Dec_1lag[[3]],Dec_2lag[[3]]), ybreaks = waiver(),
                       area_sel='to_northsea', var="Length", varlabel="Body length (cm)",
                       fix_date = TRUE, mean_date = 150)
pp3_NS <- create_panel(model=list(Dec_0lag[[3]],Dec_1lag[[3]],Dec_2lag[[3]]), ybreaks = waiver(), ylabels = waiver(),
                       area_sel='to_northsea', var="julian_recapture_std_scaled", varlabel="Recapture date (std)",
                       fix_date = TRUE, mean_date = 150)
pp4_NS <- create_panel(model=list(Dec_0lag[[3]],Dec_1lag[[3]],Dec_2lag[[3]]), ybreaks = waiver(), ylabels=waiver(),
                       area_sel='to_northsea', var="Release_year", varlabel="Release year",
                       fix_date = TRUE, mean_date = 150)
pp1_CE <- create_panel(model=list(Dec_0lag[[4]],Dec_1lag[[4]],Dec_2lag[[4]]), ybreaks = waiver(),
                       area_sel='to_ireland', var="Latitude", varlabel="Latitude (º)",
                       fix_date = TRUE, mean_date = 240, last=TRUE)
pp2_CE <- create_panel(model=list(Dec_0lag[[4]],Dec_1lag[[4]],Dec_2lag[[4]]), ybreaks = seq(-1,2,by=1),
                       area_sel='to_ireland', var="Length", varlabel="Body length (cm)",
                       fix_date = TRUE, mean_date = 240, last=TRUE)
pp3_CE <- create_panel(model=list(Dec_0lag[[4]],Dec_1lag[[4]],Dec_2lag[[4]]), ybreaks = waiver(), ylabels = waiver(),
                       area_sel='to_ireland', var="julian_recapture_std_scaled", varlabel="Recapture date (std)",
                       fix_date = TRUE, mean_date = 240, last=TRUE)
pp4_CE <- create_panel(model=list(Dec_0lag[[4]],Dec_1lag[[4]],Dec_2lag[[4]]), ybreaks = waiver(), ylabels= waiver(),
                       area_sel='to_ireland', var="Release_year", varlabel="Release year",
                       fix_date = TRUE, mean_date = 240, last=TRUE)

fake <- ggplot(data = data.frame(x = rnorm(90, 0,1),
                                 y=rnorm(90,0,1), col=sample(c("cycle 1","cycle 2","cycle 3"),90, replace=TRUE)),
               aes(x=x,y=y,col=col)) + geom_line(size=1.5) +
                scale_color_manual(name = "Migration cycle: ", values=c("red", "orange", "#FEE724")) + theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 15), legend.title = element_text(size = 16))

# The facet labels
pIC <- ggplot(data = data.frame(label="P(ICW)", x=0, y=0), aes(x=x,y=y,label=label)) + geom_text(angle =-90, size=6) + theme_void()
pNO <- ggplot(data = data.frame(label="P(NWS)", x=0, y=0), aes(x=x,y=y,label=label)) + geom_text(angle =-90, size=6) + theme_void()
pNS <- ggplot(data = data.frame(label="P(GNS)", x=0, y=0), aes(x=x,y=y,label=label)) + geom_text(angle =-90, size=6) + theme_void()
pCE <- ggplot(data = data.frame(label="P(CES)", x=0, y=0), aes(x=x,y=y,label=label)) + geom_text(angle =-90, size=6) + theme_void()


# Making the final plot
ppp <- grid.arrange(
  grobs = list(pp2,pp1,pp3,pp4,pIC,
               pp2_NO,pp1_NO,pp3_NO,pp4_NO,pNO,
               pp2_NS,pp1_NS,pp3_NS,pp4_NS,pNS,
               pp2_CE,pp1_CE,pp3_CE,pp4_CE,pCE,
               get_legend(fake)),
  heights = c(4, 4, 4, 4, 1),
  widths= c(5,5,5,5,0.7),
  layout_matrix = rbind(c( 1, 2, 3, 4, 5),
                        c( 6, 7, 8, 9,10),
                        c(11,12,13,14,15),
                        c(16,17,18,19,20),
                        c(21, 21,21,21,21)),
  left= text_grob(label="Marginal effect (scaled)", size=17, rot=90, face = "bold")
  )

ggsave(ppp, file=paste0(getwd(), "/MS/figs/Fig5.pdf"), width=16, height=10, dpi=400)











##' Function to extract the data
##' @param month the month up to which you want to keep the data for
##' @param data_origin  the migration cycle from which to extract the data
##' @details extract some limited data information for doing the boxplot in figure 4
##' @return a ggplot object to produce fig 4
##' @export create_data
##'
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
All_lags_full <- rbind(All_lags, All_lags_new) %>% mutate(direction_short = as.factor(direction)) %>%
  filter(direction_short != "All") %>%
  mutate(direction_short = recode_factor(direction_short, "Icelandic Waters" = "ICW", "Norwegian Sea"="NWS", "Greater North Sea"="GNS", "Celtic Seas"="CES"),
         cycle = as.factor(cycle),
         Release_year = as.factor(Release_year)) %>% mutate(direction_short = droplevels(direction_short))

All_lags_full <- All_lags_full %>% mutate(cycle = paste("Migration ", cycle))


nsamp <- All_lags_full %>% group_by(cycle,direction_short) %>%
  summarize(n = n(), mean = mean(Latitude), Latitude = 62) %>%
  mutate(nsamp = paste0("n=", n))

p1 <- All_lags_full %>% group_by(direction_short,cycle) %>% mutate(mediany = median(Latitude)) %>%
  ggplot(., aes(y= Latitude, x = direction_short)) +
  geom_boxplot(aes(y= Latitude, x = direction_short), outlier.color=grey(0.8), outlier.shape = 2) +
  geom_text(data=nsamp, aes(y= Latitude, x = direction_short, label=nsamp), angle=0, size=4) +
  scale_y_continuous(breaks=seq(52,64,by=4)) +
  facet_grid(. ~ cycle, scales="free_y") + xlab("ICES ecoregion") +
  theme(axis.text.x = element_text(angle = 90, size=11)) +
  ylab("Latitude (º)") + theme_bw()+
theme(axis.title = element_text(size=15), axis.text = element_text(size=12), strip.text.x = element_text(size = 15),
      axis.title.x = element_blank(), strip.background = element_rect(fill="white", color="white"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linetype = "dashed"),
      panel.grid.minor.y = element_blank())

dat <- ggplot_build(p1)$data[[1]] %>% mutate(cycle = nsamp$cycle, mean = nsamp$mean)
p1 <- p1 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=mean, yend=mean), colour="red", size=1)


nsamp1 <- All_lags_full %>% group_by(cycle,direction_short) %>%
  summarize(n = n(), mean = mean(Length), Length = 42) %>% mutate(nsamp = paste0("n=", n))

p2 <- All_lags_full %>% group_by(direction_short,cycle) %>% mutate(lengthy = median(Length)) %>%
  ggplot(.) +
  geom_boxplot(aes(y= Length, x = direction_short), outlier.color=grey(0.8), outlier.shape = 2) +
  geom_text(data=nsamp1, aes(y= Length, x = direction_short, label=nsamp), angle=0, size=4) +
  facet_grid(.  ~ cycle, scales="free_y") + xlab("ICES ecoregion") +
  theme(axis.text.x = element_text(angle = 90, size=11)) + ylab("Length (cm)")  +
  theme_bw() + coord_cartesian(ylim=c(29,42)) +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=12),
        strip.text = element_blank(),  panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.grid.minor.y = element_blank())

dat <- ggplot_build(p2)$data[[1]] %>% mutate(cycle = nsamp1$cycle, mean = nsamp1$mean)
p2 <- p2 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=mean, yend=mean), colour="red", size=1)


pp <- grid.arrange(p1,p2,nrow=2)

ggsave(pp, file="MS/Figs/Length_distribution_boxplot_summary.pdf", width=8, height=10, dpi=400)








#
# # -----------------------------------------------------------------------
# # ---- old Fig 1
# # -----------------------------------------------------------------------
#
# library(ggOceanMaps)
# library(sf)
# library(ncdf4) # package for netcdf manipulation
#
# Norway <- st_read("shapefile/ne_10m_land.shp")
# Norway <- st_crop(Norway, c(xmin = -34, ymin = 49, xmax = 22, ymax = 72))
#
# dt <- data.frame(lon = c(-31, -31, 19, 19), lat = c(45, 70, 70, 45))
#
# Catch_data <- taggart::tg_catches(lower=T)
# Mark_recap_data <- taggart::tg_expeditions(lower=T)
# d <- Mark_recap_data %>% filter(lubridate::year(releasedate) <2021, lubridate::year(releasedate) >= 2014)
# d <- filter(d,
#             ices_rectangle %in% unique(d$ices_rectangle)[unlist(sapply(31:52, function(x) grep(x, unique(d$ices_rectangle))))])
#
# d <- filter(d,
#             !(ices_rectangle %in% unique(d$ices_rectangle)[unlist(sapply("F", function(x) grep(x, unique(d$ices_rectangle))))]))
#
#
# # Releases:
# df_rel <- d %>% filter(between(lubridate::year(releasedate),2014,2020))
# df_rel <- mutate(df_rel,
#                  iceslon =mapplots::ices.rect(ices_rectangle)$lon,
#                  iceslat = mapplots::ices.rect(ices_rectangle)$lat)
#
#
# # Recaptures:
# recap_raw <- df_rel %>% filter(!is.na(recapturedate)) %>%
#   left_join(Catch_data, by = "catchid") %>%
#   mutate(RecaptureYear = ifelse(lubridate::month(recapturedate)>2,
#                                 lubridate::year(recapturedate),
#                                 lubridate::year(recapturedate)-1),
#          ReleaseYear = lubridate::year(releasedate),
#          YearsOut = RecaptureYear - ReleaseYear) %>%
#   filter(YearsOut<3, between(ReleaseYear,2014,2020),
#          !is.na(ices_rectangle.y))
#
#
# nc_data <- nc_open("data/current.nc")
# nc_data <- nc_open("data/current_new.nc")
#
# print(nc_data)
# lon <- ncvar_get(nc_data, "longitude")
# length(lon)
# lat <- ncvar_get(nc_data, "latitude")
# length(lat)
# time <- ncvar_get(nc_data, "time")
# uo <- ncvar_get(nc_data, "uo")
# dim(uo)
#
# vo <- ncvar_get(nc_data, "vo")
# dim(vo)
#
# vo_df <- reshape2::melt(vo, varnames = c("lon", "lat", "time"), value.name = "NS_comp")
# vo_df$lon <- as.vector(lon)[match(vo_df$lon, unique(vo_df$lon))]
# vo_df$lat <- as.vector(lat)[match(vo_df$lat, unique(vo_df$lat))]
# vo_df$time <- as.vector(time)[match(vo_df$time, unique(vo_df$time))]
# uo_df <- reshape2::melt(uo, varnames = c("lon", "lat", "time"), value.name = "EW_comp")
# uo_df$lon <- as.vector(lon)[match(uo_df$lon, unique(uo_df$lon))]
# uo_df$lat <- as.vector(lat)[match(uo_df$lat, unique(uo_df$lat))]
# uo_df$time <- as.vector(time)[match(uo_df$time, unique(uo_df$time))]
#
# Dat <- vo_df
# Dat$EW_comp <- uo_df$EW_comp
# Dat$date <- as.POSIXct(Dat$time*3600, origin='1950-01-01 00:00:00', tz="GMT", format = "%Y-%m-%d")
# Dat$speed <- sqrt(Dat$EW_comp^2 + Dat$NS_comp^2)
#
#
# ### May current data
# Current_data <- Dat %>% filter(time %in% c(564276, 573036, 581820, 590580, 599340, 608100, 616884))
# Current_data$Release_year <- (2014:2020)[match(Current_data$time, c(564276, 573036, 581820, 590580, 599340, 608100, 616884))]
# Current_data <- Current_data %>% filter(lat>51, lon<15, lat<68.5)
#
# survey <- Current_data %>% st_as_sf(crs = 4326, coords = c("lon", "lat"))
# surv_utm_coords <- st_coordinates(survey)
# Current_data$X <- surv_utm_coords[,1]
# Current_data$Y <- surv_utm_coords[,2]
#
# X <- basemap_data(limits = NULL, data = dt, shapefiles = NULL,
#                   bathymetry = TRUE, glaciers = FALSE, resolution = "low",
#                   lon.interval = NULL, lat.interval = NULL,
#                   rotate = FALSE, verbose = FALSE)
#
# Bathy <- X$shapefiles$bathy %>% st_as_sf(crs = 4326, coords = c("lon", "lat"))
#
# cc <- scales::seq_gradient_pal("lightblue1", "darkblue", "Lab")(seq(0,1,length.out=8))
#
# new_sf <- ICESecoregion %>% filter(Ecoregion %in% c("Norwegian Sea", "Icelandic Waters", "Greater North Sea", "Celtic Seas", "Faroes")) %>%
#   group_by(Ecoregion) %>% sfheaders::sf_remove_holes()
#
#
#
# g1 <- ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
#   geom_sf(data = Bathy, aes(col=depth, fill=depth)) +
#   scale_colour_manual(values=cc, name = "Depth") + scale_fill_manual(values=cc, name = "Depth") +
#   geom_sf(data=Norway, col=grey(0.3)) +
#   # new_scale_color() +
#   geom_sf(data=new_sf, aes(geometry = geometry), col = "black", fill=NA, size=0.9) +
#   geom_text(data = data.frame(x=c(-19,0,3.2,-13.3,-6.8), y=c(61.3,69.5,55.8,55,62.5),
#                               label=c("Icelandic \nWaters", "Norwegian Sea", "Greater \nNorth \nSea", "Celtic\nSeas", "Faroes")),
#             aes(x=x, y=y, label=label), col="white", size=5, fontface=2) +
#   geom_point(data = df_rel %>% distinct(longitude,latitude),    aes(x = longitude, y = latitude), col = "yellow", alpha = 0.5, size = .5)+
#   geom_point(data = recap_raw %>% distinct(clon,clat), aes(x = clon, y = clat), col = "red", size = .6)+
#   coord_sf(xlim=c(-31, 19), ylim=c(50, 70)) +
#   labs(x="Longitude", y="Latitude") +
#   theme_bw() + theme(plot.title = element_text(hjust = 0.5),
#                      legend.title = element_text(face = "bold", size = 12),
#                      strip.background = element_rect(fill = "white"),
#                      strip.text = element_text(face ="bold", size = 14),
#                      axis.text = element_text(size =13),
#                      axis.title = element_text(size =14))
# #panel.spacing.x = unit(6.5, "mm"))+
# # ggsave(g1, filename = "MS/figs/Fig2_May_Current_bathym.pdf",
# #        width=28, height=32, units="cm", dpi = 450)
# #
#
# #
# #   	ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
# #     	  geom_sf(data=bla, aes(geometry = geometry, col=Ecoregion), fill=NA, size=2)
#
#
# g2 <- ggplot(data=Norway) + geom_sf(col=grey(0.3)) +
#   geom_sf(data = Bathy, aes(col=depth, fill=depth), show.legend = FALSE) +
#   scale_colour_manual(values=cc, name = "Depth") + scale_fill_manual(values=cc, name = "Depth") + geom_sf(data=Norway, col=grey(0.3)) +
#   coord_sf(xlim=c(-31, 19), ylim=c(50, 70)) +
#   new_scale_color() +
#   new_scale_fill() +
#   geom_segment(data = Current_data %>% filter(speed > 0.15, Release_year %in% c(2014, 2018)), aes(x=X, y=Y, xend = X + EW_comp*5, yend = Y + NS_comp*5, col=speed),
#                arrow = arrow(angle = 15, length = unit(0.03, "inches"), type = "closed"), alpha = 0.8, size=0.3) +
#   scale_color_gradient(low="lightpink1", high="darkred", name = "Current\nspeed (m/s)") +
#   facet_wrap(~ Release_year, ncol=2) + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
#                                                           legend.title = element_text(face = "bold", size = 12),
#                                                           strip.background = element_rect(fill = "white"),
#                                                           strip.text = element_text(face ="bold", size = 14),
#                                                           axis.text = element_text(size =13),
#                                                           axis.title = element_text(size =14)) +
#   labs(x="Longitude", y="Latitude")
#
#
# pp <- gridExtra::grid.arrange(g1, g2, nrow=2)
# ggsave(pp, file="MS/Figs/Fig1.png", width=10, height=10, dpi=400)
#
#
#
#
#


