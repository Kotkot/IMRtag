N_cross = 500
P_cross <- 0.8

# Just in case you have not compiled and loaded the model
# use_version_simple <- paste0(getwd(), "/src/mackerel_mvt_model_nothresh")
# compile(paste0(use_version_simple, ".cpp"))
# dyn.load(use_version_simple)
#
# use_version <- paste0(getwd(), "/src/mackerel_mvt_model")
# compile(paste0(use_version, ".cpp"))
# dyn.load(use_version)
#

cross_validation <- function(x){
  set.seed(x)

      N_threshold = 1
      training <- Data_mackerel_use_Ireland_select %>% group_by(Release_year) %>% slice_sample(prop=P_cross) %>% ungroup()
      testing <- Data_mackerel_use_Ireland_select[-training$ID,]

      m1 <- lm(log_rate ~ factor(Tag_area)*factor(Release_timing) + factor(Release_year) + length_scaled + julian_recapture_scaled, data=training)
      m_frame <- model.frame(m1)
      XX <- model.matrix(m1, m_frame)
      XX_pred <- model.matrix(m1, data=testing)

    # Running the no-threshold model
      data_tmb <- list(N=nrow(training),   # number of data points
                       X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
                       Nthres=length(threshold_vals),
                       y = training$log_rate,
                       X_pred=as.matrix(as.data.frame(XX_pred)),
                       N_pred=nrow(testing),
                       Likconfig = 0)      # 0 = dnorm, 1 = dgamma
      parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                           log_sigma = rep(log(0.2),N_threshold))
      Map = list()
      obj <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model_nothresh", map=Map)
      opt <- fit_tmb( obj=obj, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000), quiet=TRUE)
      sd_report <- sdreport(obj)
      mu_pred <- matrix(summary(sd_report, "report")[grep("mu_pred", rownames(summary(sd_report, "report"))),1], ncol=1, byrow=FALSE)

      mse_nothres <- mean((mu_pred - testing$log_rate)^2)

    # Running the best threshold model
      N_threshold <- 2
      data_tmb1 <- list(K=N_threshold,  # number of mixture components
                       N=nrow(training),   # number of data points
                       X=as.matrix(as.data.frame(XX)),          # the design matrix for the fixed effect
                       Nthres=length(threshold_vals),
                       thresh=threshold_vals,
                       mean_diff_tag_area= mean_diff_tag_area,
                       is_from_west=ifelse(training$Tag_area == "West_Ireland",1,0),
                       thres_cov = training$julian_std,
                       y = training$log_rate,
                       X_pred=as.matrix(as.data.frame(XX_pred)),
                       N_pred=nrow(testing),
                       Likconfig = 0)      # 0 = dnorm, 1 = dgamma

      best_mod <- c(which.min(AIC_selection),which.min(AIC_selection2))[which.min(c(AIC_selection[which.min(AIC_selection)],AIC_selection2[which.min(AIC_selection2)]))]

      Map1 = list()
      Map1$beta <- factor(Maps_best[best_mod,])
      Map1$log_sigma <- factor(c(2*ncol(XX)+1:2))
      parameters_tmb1 <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                             log_sigma = rep(log(0.2),N_threshold))
      parameters_tmb1$beta[which(is.na(Maps_best[best_mod,]))] <- 0

      obj1 <- MakeADFun(data_tmb1, parameters_tmb1, random = NULL, DLL = "mackerel_mvt_model", map=Map1)
      opt1 <- fit_tmb( obj=obj1, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000), quiet=TRUE)
      sd_report1 <- sdreport(obj1)
      mu_pred <- matrix(summary(sd_report1, "report")[grep("mu_pred", rownames(summary(sd_report1, "report"))),1], ncol=2, byrow=FALSE)

      LL <- obj1$report()$LL

      weight <- apply(exp(LL), 2, mean)
      weight <- weight/sum(weight)

      Prediction <- rep(0, data_tmb1$N_pred)
      for (n in 1:data_tmb1$N_pred){
        for (thr in 1:data_tmb1$Nthres){
          if (testing$julian_std[n] < (data_tmb1$thresh[thr])){
            Prediction[n] = Prediction[n] + weight[thr]*mu_pred[n,1]
          }
          if (testing$julian_std[n] >= (data_tmb1$thresh[thr])){
            Prediction[n] = Prediction[n] + weight[thr]*mu_pred[n,2]
          }
        }
      }
      mse_thres <- mean((Prediction - testing$log_rate)^2)

    # Output the MSE of the model types
      out <- c(mse_nothres, mse_thres)
      return(out)
 }

Cross_val_res <- 1:N_cross %>% map(function(x) cross_validation(x))
Cross_val_res_df <- do.call(rbind, Cross_val_res)
colnames(Cross_val_res_df) <- c("No_thres", "Thresh")
boxplot(Cross_val_res_df)

Cross_val_res_df_melt <- reshape2::melt(Cross_val_res_df)
ggplot(Cross_val_res_df_melt, aes(x=Var2, y=value)) + geom_violin() + geom_boxplot()

t.test(x=Cross_val_res_df[,1], y=Cross_val_res_df[,2], alternative = "two.sided")
