### Testing the random effect model

  ## the no-trehsold model
    use_version_simple <- paste0(getwd(), "/src/mackerel_mvt_model_nothresh_RE")
    compile(paste0(use_version_simple, ".cpp"))
    dyn.load(use_version_simple)

    N_threshold <- 1
    data_tmb <- list(N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
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

    parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                           log_sigma = log(0.2),
                           year = rep(0, length(unique(Data_mackerel_use_Ireland_select$Release_year))),
                           log_sigma_year = log(0.2)
    )

    Map = list(beta = 1:(ncol(XX)),
               log_sigma = ncol(XX)+1,
               year=ncol(XX)+1+1:length(parameters_tmb$year),
               log_sigma_year = ncol(XX)+8)
    Map$beta <- factor(Map$beta)
    Map$log_sigma <- factor(Map$log_sigma)
    Map$year <- factor(Map$year)
    Map$log_sigma_year <- factor(Map$log_sigma_year)

    obj <- MakeADFun(data_tmb, parameters_tmb, random = "year", DLL = "mackerel_mvt_model_nothresh_RE", map=Map)
    opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

    opt


  ## The threshold model
    use_version <- paste0(getwd(), "/src/mackerel_mvt_model_RE")
    compile(paste0(use_version, ".cpp"))
    dyn.load(use_version)

    m1 <- lm(log_rate ~ factor(Tag_area)*factor(Release_timing) + length_scaled + julian_recapture_scaled, data=Data_mackerel_use_Ireland_select)
    m_frame <- model.frame(m1)
    XX <- model.matrix(m1, m_frame)


    N_threshold <- 2
    data_tmb <- list(K=N_threshold,  # number of mixture components
                     N=nrow(Data_mackerel_use_Ireland_select),   # number of data points
                     X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
                     Nthres=length(threshold_vals),
                     thresh=threshold_vals,
                     mean_diff_tag_area= mean_diff_tag_area,
                     is_from_west=ifelse(Data_mackerel_use_Ireland_select$Tag_area == "West_Ireland",1,0),
                     thres_cov = Data_mackerel_use_Ireland_select$julian_std,
                     y = Data_mackerel_use_Ireland_select$log_rate,
                     X_pred =as.matrix(as.data.frame(XX)),
                     N_pred =nrow(Data_mackerel_use_Ireland_select),
                     N_year = length(unique(Data_mackerel_use_Ireland_select$Release_year)),
                     Year_ID = as.numeric(as.character(Data_mackerel_use_Ireland_select$Release_year))-2014,
                     Year_ID_pred = as.numeric(as.character(Data_mackerel_use_Ireland_select$Release_year))-2014,
                     Likconfig = 0      # 0 = dnorm, 1 = dgamma
    )

    parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                           log_sigma = rep(log(0.2),N_threshold),
                           year = matrix(rep(0, length(unique(Data_mackerel_use_Ireland_select$Release_year))*N_threshold), ncol=2),
                           log_sigma_year = rep(log(0.2),N_threshold)
    )

    Map = list(beta = matrix(1:(ncol(XX)*2), ncol = 2),
               log_sigma = c(2*ncol(XX)+1:2),
               log_sigma_year = c(2*ncol(XX)+3:4))
    Map$beta <- factor(Map$beta)
    Map$log_sigma <- factor(Map$log_sigma)
    Map$log_sigma_year <- factor(Map$log_sigma_year)

    obj1break <- MakeADFun(data_tmb, parameters_tmb, random = c("year"), DLL = "mackerel_mvt_model_RE", map=Map)
    opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

    sde <- sdreport(obj1break)
    summary(sde, "fixed")
    summary(sde, "random")
    mu_pred <- matrix(summary(sde, "report")[grep("mu_pred", rownames(summary(sde, "report"))),1], ncol=2, byrow=FALSE)

    # Calculating the actual prediction
    LL <- obj1break$report()$LL

    # The weighting factor is the average likelihood across the N datapoint per threshold value
    weight <- apply(exp(LL), 2, mean)
    weight <- weight/sum(weight)
    plot(data_tmb$thresh, weight, type="b")
    #           weight <- rep(0,17)
    # 					weight[5]=1

    Prediction <- rep(0, data_tmb$N)
    for (n in 1:data_tmb$N){
      for (thr in 1:data_tmb$Nthres){
        if (data_tmb$thres_cov[n] < (data_tmb$thresh[thr])){
          Prediction[n] = Prediction[n] + weight[thr]*mu_pred[n,1]
        }
        if (data_tmb$thres_cov[n] >= (data_tmb$thresh[thr])){
          Prediction[n] = Prediction[n] + weight[thr]*mu_pred[n,2]
        }
      }
    }

    par(mfrow=c(2,1), mar=c(4,3,1,1), oma=c(1,1,1,1))
    plot(Prediction, data_tmb$y); abline(0,1)
    qqnorm(y=(Prediction-data_tmb$y)/sd(Prediction-data_tmb$y), xlim=c(-3.5,3.5),ylim=c(-3.5,3.5))
    abline(0,1, lty=2)

    opt1break$objective * 2 + p * length(opt1break$par)
