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
  library(DHARMa)
  library(ggeffects)
  library(groupdata2)
  library(furrr)
  tonum <- function(x) as.numeric(as.character(x))

#### loadig some shapefiles
  new_proj <- 3035
  Norway <- st_read("D:/Dropbox/IMR_projects/Shapefiles/ne_10m_land.shp")
  Norway <- st_crop(Norway, c(xmin = -35, ymin = 51, xmax = 21, ymax = 73))
  Norway_proj <- st_transform(Norway, crs=new_proj)
  ICESecoregion <- st_read("D:/Dropbox/IMR_projects/Shapefiles/ICES_ecoregions_20171207_erase_ESRI.shp")
  ICESecoregion_cut <- st_crop(ICESecoregion, c(xmin = -35, ymin = 51, xmax = 21, ymax = 73))
  ICESecoregion_proj <- st_transform(ICESecoregion, crs=new_proj)
  ICESecoregion_cut_proj <- st_transform(ICESecoregion_cut, crs=new_proj)

  ggplot(ICESecoregion_cut_proj) + geom_sf(aes(geometry = geometry, fill=Ecoregion))

#### Now load the data for the analysis
  cutoff_months = c(6:14)
  source("R/0_prepare_data.R")

#### Some plotting of the data
  ggplot(Norway) + geom_sf() + geom_jitter(data=Data_mackerel_use_Ireland_select_origin, aes(x=cLon, y=cLat, col=factor(Catch_month))) +
    theme_bw() + facet_wrap(~Catch_year) + geom_hline(yintercept=62) + geom_vline(xintercept=-10)+ geom_vline(xintercept=-4, col="red")
  ggplot(Norway) + geom_sf() + geom_jitter(data=Data_mackerel_use_Ireland_select_origin, aes(x=Longitude, y=Latitude, col=factor(Catch_month))) +
    theme_bw() + facet_wrap(~Catch_year) + geom_hline(yintercept=62) + geom_vline(xintercept=-10)+ geom_vline(xintercept=-4, col="red") +
    coord_sf(ylim=c(50, 60), xlim=c(-15,0))

#### Now load the functions for the analysis
  source("R/Functions.R")

#### Then run the analysis
    Dec_0lag <- run_directionality(month=14, data_origin="cycle1", model_selection = "none", scale = "linear", alpha = 0.05)
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle1_IS.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_0lag[[1]][[9]], dat=Dec_0lag[[1]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle1_NO.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_0lag[[2]][[9]], dat=Dec_0lag[[2]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle1_NS.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_0lag[[3]][[9]], dat=Dec_0lag[[3]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle1_IR.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_0lag[[4]][[9]], dat=Dec_0lag[[4]]$data)
  	dev.off()

  	Dec_1lag <- run_directionality(month=14, data_origin="cycle2", model_selection = "none", scale = "linear")
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle2_IS.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_1lag[[1]][[9]], dat=Dec_1lag[[1]]$data, nsim=499)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle2_NO.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_1lag[[2]][[9]], dat=Dec_1lag[[2]]$data, nsim=499)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle2_NS.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_1lag[[3]][[9]], dat=Dec_1lag[[3]]$data, nsim=5000)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle2_IR.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_1lag[[4]][[9]], dat=Dec_1lag[[4]]$data, nsim=5000)
  	dev.off()

  	Dec_2lag <- run_directionality(month=14, data_origin="cycle3", model_selection = "none",scale = "linear")
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle3_IS.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_2lag[[1]][[9]], dat=Dec_2lag[[1]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle3_NO.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_2lag[[2]][[9]], dat=Dec_2lag[[2]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle3_NS.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_2lag[[3]][[9]], dat=Dec_2lag[[3]]$data)
  	dev.off()
  	png(filename=paste0(getwd(), "/MS/figs/Resid_cycle3_IR.png"), res=400, width=16, height=12, units="cm")
  	Resid_plot(mod=Dec_2lag[[4]][[9]], dat=Dec_2lag[[4]]$data)
  	dev.off()

  # main models based on model selection
  	Dec_0lag_AIC <- run_directionality(month=14, data_origin="cycle1", model_selection = "AIC",  scale = "linear")
  	Dec_1lag_AIC <- run_directionality(month=14, data_origin="cycle2", model_selection = "AIC",  scale = "linear")
  	Dec_2lag_AIC <- run_directionality(month=14, data_origin="cycle3", model_selection = "AIC",  scale = "linear")


	# Sensitivity
  	Jan_0lagb <- run_directionality(month=13, data_origin="cycle1", model_selection = "none", scale = "linear")
  	Jan_1lagb <- run_directionality(month=13, data_origin="cycle2", model_selection = "none", scale = "linear")
  	Jan_2lagb <- run_directionality(month=13, data_origin="cycle3", model_selection = "none", scale = "linear")



