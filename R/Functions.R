##' Function to create the marginal effect plot
##' @param it the iteratiin number (this sets the seed)
##' @param kfolds the number of k-fold to perform
##' @param model_type  the "go-to-ecoregion" to focus the analysis on
##' @param inputdata  the corresponding data
##' @details
##' @return a list containing the fitted model objects, their AIC, BIC, and MSE values, and the data
##' @export kfolding
##'
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
    m02 <- (gam(y ~  -1 + Release_year+ s(Latitude,k=3), family=binomial,
                data = dat_use))
    m03 <- (gam(y ~  -1 + Release_year+ s(Length), family=binomial,
                data = dat_use))
    m04 <- (gam(y ~  -1 + Release_year+ s(julian_recapture_std_scaled, k=3), family=binomial,
                data = dat_use))
    m05 <- (gam(y ~  -1 + Release_year+ s(Latitude,k=3) + s(julian_recapture_std_scaled, k=3), family=binomial,
                data = dat_use))
    m06 <- (gam(y ~  -1 + Release_year+ s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
                data = dat_use))
    m07 <- (gam(y ~  -1 + Release_year+ s(Latitude,k=3) + s(Length), family=binomial,
                data = dat_use))
    if (month <= 12) m1 <- (gam(y ~  -1 + s(Latitude,k=3) + Release_year  + s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
                                data = dat_use))
    if (month > 12) m1 <- (gam(y ~  -1 + s(Latitude,k=3) + Release_year  + s(Length,k=3) + s(julian_recapture_std_scaled, k=3), family=binomial,
                               data = dat_use))
    # m2 <- (gam(y ~  -1 + s(Latitude) + Release_year + s(Length, by=Release_year) + s(julian_recapture_std_scaled, k=3), family=binomial,
    #            data = dat_use))
    # m3 <- (gam(y ~  -1 + s(Latitude, by=Release_year) + Release_year + s(Length) + s(julian_recapture_std_scaled, k=3), family=binomial,
    #            data = dat_use))
    # m4 <- (gam(y ~  -1 + s(Latitude) + Release_year  + s(Length) + s(julian_recapture_std_scaled, by=Release_year), family=binomial,
    #            data = dat_use))
    m2 <- (gam(y ~  -1 + s(Latitude,k=3) + Release_year + s(julian_recapture_std_scaled, by=length_bin), family=binomial,
               data = dat_use))
    m3 <- (gam(y ~  -1 + s(Latitude,k=3) + Release_year + s(Length, k=3) + s(julian_recapture_std_scaled, by=length_bin), family=binomial,
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


##' Function to create the marginal effect plot
##' @param month the month to keep the data up-to
##' @param data_origin is the "go-to" direction to choose: choice among the four ICES ecoregions
##' @param model_selection  whether to perform a model selection "none" or "AIC"
##' @param scale  the labels for the y axis. needed for plotting size adjustement
##' @param alpha  to you want to fix the date for the marginal effect plot? if so specify in mean_date the date to use
##' @details creates all the individual panel that aggregates results across migration cycles
##' @return a list of model runs for all ecoregions + produced figures
##' @export run_directionality
##'
run_directionality <- function(month, data_origin="cycle1", model_selection = "none",
                               sensitivity=FALSE, scale = "response", alpha= 0.05, adjust_ecoregion = FALSE){
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

  intersection <- st_intersects(Data_mackerel_use_Ireland_select_sf, ICESecoregion_proj, sparse=FALSE)
  area <- rep(NA, nrow(intersection))
  for (i in 1:nrow(intersection)) { if (length(which(intersection[i,]==TRUE)>0)) area[i]= which(intersection[i,]==TRUE) }
  # a little adjustment of the area falling on land
  area[which(is.na(area == TRUE) & (Data_mackerel_use_Ireland_select[, 'cLon'] == -1.5))] = 9
  area[which(is.na(area == TRUE) & (Data_mackerel_use_Ireland_select[, 'cLon'] == 5.5))] = 11
  # points falling in faroes but re-allocated to norway and iceland
  area[which(area == 15 & Data_mackerel_use_Ireland_select$cLon > -5)] = 16
  area[which(area == 15 & Data_mackerel_use_Ireland_select$cLon < -5)] = 13
  Data_mackerel_use_Ireland_select$direction <- ICESecoregion_proj$Ecoregion[area]
  Data_mackerel_use_Ireland_select$ID <- 1:nrow(Data_mackerel_use_Ireland_select)

  # Rescaling parameters to ease interpretation
  Data_mackerel_use_Ireland_select$julian_recapture_scaled <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std)
  Data_mackerel_use_Ireland_select$julian_recapture_standardized <- scale(Data_mackerel_use_Ireland_select$julian_recapture_std, center=FALSE)
  Data_mackerel_use_Ireland_select$length_scaled <- scale(Data_mackerel_use_Ireland_select$Length)
  Data_mackerel_use_Ireland_select$Latitude_scaled <- scale(Data_mackerel_use_Ireland_select$Latitude)
  Data_mackerel_use_Ireland_select$Latitude2_scaled <- scale((Data_mackerel_use_Ireland_select$Latitude)^2)

  ##Adjustment of the directionality
  if (adjust_ecoregion == TRUE) {
    Data_mackerel_use_Ireland_select <- Data_mackerel_use_Ireland_select %>%
    mutate(direction= ifelse((direction == "Celtic Seas" & cLat>62), "Norwegian Sea", direction))
  }


  ## Creating binary response variable
  Data_mackerel_use_Ireland_select <- Data_mackerel_use_Ireland_select %>%
    mutate(toiceland= ifelse(direction == "Icelandic Waters", 1, 0),
           to_norway = ifelse(direction == "Norwegian Sea", 1, 0),
           to_northsea = ifelse(direction == "Greater North Sea", 1, 0),
           to_ireland = ifelse(direction == "Celtic Seas", 1, 0)
    )



  ## Analysis of the movement directionality
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

    ### Go to Iceland
    main_iceland = kfolding(it=1, kfolds=1, model_type="Iceland", inputdata=www)
    if (model_selection == "AIC") best_iceland <- main_iceland[[which.min(main_iceland$AICs$AIC[1:9])]]
    if (model_selection == "none") best_iceland <- main_iceland[[9]]
    simul <- simulateResiduals(best_iceland, plot=TRUE)
    par(mfrow=c(2,2))
    plot(www$Latitude, simul$scaledResiduals)
    plot(www$Length, simul$scaledResiduals)
    plot(www$julian_recapture_std, simul$scaledResiduals)
    plot(www$Release_year, simul$scaledResiduals)


    ### Go to norway
    main_norway = kfolding(it=1, kfolds=1, model_type="Norway", inputdata=www)
    if (model_selection == "AIC") best_norway <- main_norway[[which.min(main_norway$AICs$AIC[1:9])]]
    if (model_selection == "none") best_norway <- main_norway[[9]]
    simul <- simulateResiduals(best_norway, plot=TRUE)
    par(mfrow=c(2,2))
    plot(www$Latitude, simul$scaledResiduals)
    plot(www$Length, simul$scaledResiduals)
    plot(www$julian_recapture_std, simul$scaledResiduals)
    plot(www$Release_year, simul$scaledResiduals)

    ### Go to north sea
    main_northsea = kfolding(it=1, kfolds=1, model_type="Northsea", inputdata=www)
    if (model_selection == "AIC") best_northsea <- main_northsea[[which.min(main_northsea$AICs$AIC[1:9])]]
    if (model_selection == "none") best_northsea <- main_northsea[[9]]
    simul <- simulateResiduals(best_northsea, plot=TRUE)
    par(mfrow=c(2,2))
    plot(www$Latitude, simul$scaledResiduals)
    plot(www$Length, simul$scaledResiduals)
    plot(www$julian_recapture_std, simul$scaledResiduals)
    plot(www$Release_year, simul$scaledResiduals)

    ### Go to Ireland
    main_ireland = kfolding(it=1, kfolds=1, model_type="Ireland", inputdata=www)
    if (model_selection == "AIC") best_ireland <- main_ireland[[which.min(main_ireland$AICs$AIC[1:9])]]
    if (model_selection == "none") best_ireland <- main_ireland[[9]]
    simul <- simulateResiduals(best_ireland, plot=TRUE)
    par(mfrow=c(2,2))
    plot(www$Latitude, simul$scaledResiduals)
    plot(www$Length, simul$scaledResiduals)
    plot(www$julian_recapture_std, simul$scaledResiduals)
    plot(www$Release_year, simul$scaledResiduals)


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
    pp_IS <- visreg::visreg(fit=best_iceland, xvar="Latitude", plot=FALSE, data=main_iceland$data, scale = scale, alpha = alpha)
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
    pp_NO <- visreg::visreg(best_norway, "Latitude", plot=FALSE, data=main_norway$data, scale = scale, alpha = alpha)
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
    pp_NS <- visreg::visreg(best_northsea, "Latitude", plot=FALSE, data=main_northsea$data, scale = scale, alpha = alpha)
    p0_NS <- ggplot(pp_NS$fit, aes(x=Latitude, y=visregFit)) + geom_line() +
      geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
      theme_bw() +
      labs(y="", x="Latitude (º)")
    if( Limit_month > 12 & sensitivity == FALSE) p0_NS <- p0_NS + ggctr2
    pp0_NS <- p0_NS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=Latitude)) +
      geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=Latitude), col="red") #+
    #geom_text(data=nsamp_NS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  } else {
    pp0_NS <- ggplot() + theme_void()
    if ( Limit_month > 12 & sensitivity == FALSE) pp0_NS <- pp0_NS + ggctr4
  }
  if ( Limit_month > 12 & sensitivity == FALSE){
    if (length(grep("Latitude", best_ireland$call))>0) {
      pp_IR <- visreg::visreg(best_ireland, "Latitude", plot=FALSE, data=main_ireland$data, scale = scale, alpha = alpha)
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
    pp_IS <- visreg::visreg(best_iceland, "Length", plot=FALSE, data=main_iceland$data, scale = scale, alpha = alpha)
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
    pp_NO <- visreg::visreg(best_norway, "Length", plot=FALSE, data=main_norway$data, scale = scale, alpha = alpha)
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
    pp_NS <- visreg::visreg(best_northsea, "Length", plot=FALSE, data=main_northsea$data, scale = scale, alpha = alpha)
    p0_NS <- ggplot(pp_NS$fit, aes(x=Length, y=visregFit)) + geom_line() +
      geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
      theme_bw() +
      labs(y="Marginal effect", x="Length (cm)")
    if( Limit_month > 12 & sensitivity == FALSE) p0_NS <- p0_NS + ggctr0 + theme(axis.title.x = element_blank())
    pp1_NS <- p0_NS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=Length)) +
      geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=Length), col="red")
  } else {
    pp1_NS <- ggplot() + theme_void()
    if ( Limit_month > 12 & sensitivity == FALSE) pp1_NS <- pp1_NS + ggctr4
  }
  if ( Limit_month > 12 & sensitivity == FALSE){
    if (length(grep("Length", best_ireland$call))>0) {
      pp_IR <- visreg::visreg(best_ireland, "Length", plot=FALSE, data=main_ireland$data, scale = scale, alpha = alpha)
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
    pp_IS <- visreg::visreg(best_iceland, "julian_recapture_std_scaled", plot=FALSE, data=main_iceland$data, scale = scale, alpha = alpha)
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
    pp_NO <- visreg::visreg(best_norway, "julian_recapture_std_scaled", plot=FALSE, data=main_norway$data, scale = scale, alpha = alpha)
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
    pp_NS <- visreg::visreg(best_northsea, "julian_recapture_std_scaled", plot=FALSE, data=main_northsea$data, scale = scale, alpha = alpha)
    p0_NS <- ggplot(pp_NS$fit, aes(x=julian_recapture_std_scaled, y=visregFit)) + geom_line() +
      geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), fill=grey(0.5), alpha=0.3) +
      theme_bw()+ labs(y="Marginal effects", x="Recapture date std")
    # scale_y_continuous(labels=func1space)
    if( Limit_month > 12 & sensitivity == FALSE) p0_NS <- p0_NS + ggctr0 + theme(axis.title.x = element_blank())
    pp2_NS <- p0_NS + geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==0), aes(x=julian_recapture_std_scaled)) +
      geom_rug(data=data.frame(www, visregFit=-Inf) %>% filter(tonorthsea_bin==1), aes(x=julian_recapture_std_scaled), col="red")# +
    #geom_text(data=nsamp_NS, aes(x=-Inf, y=Inf, label=ntot), vjust=1.2, hjust=-0.1, size=3)
  } else {
    pp2_NS <- ggplot() + theme_void()
    if ( Limit_month > 12 & sensitivity == FALSE) pp2_NS <- pp2_NS + ggctr4
  }
  if ( Limit_month > 12 & sensitivity == FALSE){
    if (length(grep("julian_recapture_std_scaled", best_ireland$call))>0) {
      pp_IR <- visreg::visreg(best_ireland, "julian_recapture_std_scaled", plot=FALSE, data=main_ireland$data, scale = scale, alpha = alpha)
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
    pp <- visreg::visreg(best_iceland, "Release_year", plot=FALSE, data=main_iceland$data, alpha = alpha)
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
    pp <- visreg::visreg(best_norway, "Release_year", plot=FALSE, data=main_norway$data, alpha = alpha)
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
    pp <- visreg::visreg(best_northsea, "Release_year", plot=FALSE, data=main_northsea$data, alpha = alpha)
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
      pp_IR <- visreg::visreg(best_ireland, "Release_year", plot=FALSE, data=main_ireland$data, alpha = alpha)
      pp$fit <- pp$fit %>% mutate(Release_year_fct = factor(Release_year, levels=2014:2020))
      pp3_IR <- ggplot(pp$fit, aes(x=Release_year_fct, y=visregFit)) + geom_point() +
        geom_errorbar(aes(ymin=visregLwr, ymax=visregUpr), col=grey(0.5), width=0.4) +
        theme_bw() +scale_x_discrete(drop=FALSE) +
        ylab("") +
        xlab("Year") + labs(tag="P(Celtic Seas)") +
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
  ppp <- grid.arrange(pp0_IS,pp1_IS,pp2_IS,pp3_IS,
                        pp0_NO,pp1_NO,pp2_NO,pp3_NO,
                        pp0_NS,pp1_NS,pp2_NS,pp3_NS,
                        pp0_IR,pp1_IR,pp2_IR,pp3_IR,
                        layout_matrix=matrix(c(1:16), byrow=T, ncol=4, nrow=4), widths=c(1.07,1,1,1.1), heights=c(1,1,1,1.07))
    ggsave(ppp, file=paste0(getwd(), "/MS/figs/Marginal_cutoffmonth", month, "_lag", data_origin, model_selection, ".pdf"),width=26, height=20, units="cm", dpi = 450)


  best_models <- list(main_iceland, main_norway, main_northsea, main_ireland, ppp)

  return(best_models)
}


##' a slighty modified version of the ensureDHARMa function from the DHARMa package
##' @param simulationOutput is the output from "simulateResiduals" function in DHARMa
##' @param convert is the "go-to" direction to choose: choice among the four ICES ecoregions
##' @details a slighty modified version of the ensureDHARMa function from the DHARMa package (needed for the residual test plotQQunif1)
##' @return
##' @export ensureDHARMa1
##'
ensureDHARMa1 <- function(simulationOutput,
                          convert = F){

  if(inherits(simulationOutput, "DHARMa")){
    return(simulationOutput)
  } else {

    if(convert == FALSE) stop("wrong argument to function, simulationOutput must be a DHARMa object!")
    else {

      if (class(simulationOutput)[1] %in% getPossibleModels()){
        if (convert == "Model" | convert == T) return(simulateResiduals(simulationOutput))
      } else if(is.vector(simulationOutput, mode = "numeric") & convert == T) {
        out = list()
        out$scaledResiduals = simulationOutput
        out$nObs = length(out$scaledResiduals)
        class(out) = "DHARMa"
        return(out)
      }
    }
  }
  stop("wrong argument to function, simulationOutput must be a DHARMa object or a numeric vector of quantile residuals!")
}

##' a slighty modified version of the plotQQunif function from the DHARMa package
##' @param simulationOutput is the output from "simulateResiduals" function in DHARMa
##' @param testUniformity do uniformity test?
##' @param testOutliers  do outlier test?
##' @param testDispersion  do dispersion test?
##' @param main  the title of the plot
##' @details a slighty modified version of the plotQQunif function from the DHARMa package
##' @return
##' @export plotQQunif1
##'
plotQQunif1 <- 	function (simulationOutput, testUniformity = T, testOutliers = T,
                          testDispersion = T, main ="QQ plot residuals", ...)
{
  simulationOutput = ensureDHARMa1(simulationOutput, convert = "Model")
  gap::qqunif(simulationOutput$scaledResiduals, pch = 2, bty = "n",
              logscale = F, col = "black", cex = 0.6, main = main,
              cex.main = 1, ...)
  if (testUniformity == TRUE) {
    temp = testUniformity(simulationOutput, plot = F)
    legend("topleft", c(paste("KS test: p=", round(temp$p.value,
                                                   digits = 5)), paste("Deviation ", ifelse(temp$p.value <
                                                                                              0.05, "significant", "n.s."))), text.col = ifelse(temp$p.value <
                                                                                                                                                  0.05, "red", "black"), bty = "n")
  }
  if (testOutliers == TRUE) {
    temp = testOutliers(simulationOutput, plot = F)
    legend("bottomright", c(paste("Outlier test: p=", round(temp$p.value,
                                                            digits = 5)), paste("Deviation ", ifelse(temp$p.value <
                                                                                                       0.05, "significant", "n.s."))), text.col = ifelse(temp$p.value <
                                                                                                                                                           0.05, "red", "black"), bty = "n")
  }
  if (testDispersion == TRUE) {
    temp = testDispersion(simulationOutput, plot = F)
    legend("center", c(paste("Dispersion test: p=", round(temp$p.value,
                                                          digits = 5)), paste("Deviation ", ifelse(temp$p.value <
                                                                                                     0.05, "significant", "n.s."))), text.col = ifelse(temp$p.value <
                                                                                                                                                         0.05, "red", "black"), bty = "n")
  }
}

##' Function to create the marginal effect plot
##' @param mod is model run
##' @param dat is the corresponding data
##' @param nsim  numbre of simulation to perform (input for the simulateResiduals function in DHARMa)
##' @details creates residual QQplot as in DHARMa
##' @return creates residual QQplot as in DHARMa
##' @export Resid_plot
##'
Resid_plot <- function(mod, dat, nsim=5000){
  simul <- simulateResiduals(mod, plot=FALSE, n=nsim)
  par(mfrow=c(2,3), mar=c(4,4,1,1), oma=c(0,1,1,0))
  plotQQunif1(simul, main = "")
  dat$residuals <- simul$scaledResiduals
  dat$Release_year <- droplevels(dat$Release_year)
  plot(dat$Length, dat$residuals, xlab="Fish bodysize (cm)", ylab="Scaled residuals");# lines(dat$Length, predict(gam(residuals ~ s(Length), data=dat)))
  plot(dat$julian_recapture_std, dat$residuals, xlab="Recapture date", ylab="Scaled residuals"); #lines(dat$julian_recapture_std, predict(gam(residuals ~ s(julian_recapture_std), data=dat)))
  plot(dat$Release_year, dat$residuals, xlab="Release year", ylab="Scaled residuals");
  plot(dat$julian_release_std, dat$residuals, xlab="Release date", ylab="Scaled residuals"); #lines(dat$julian_release_std, predict(gam(residuals ~ s(julian_release_std), data=dat)))
}
