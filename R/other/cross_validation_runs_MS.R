N_cross = 500
P_cross <- 0.8
do_parallel <- TRUE

# Just in case you have not compiled and loaded the model
# use_version_simple <- paste0(getwd(), "/src/mackerel_mvt_model_nothresh_RE")
# compile(paste0(use_version_simple, ".cpp"))
# dyn.load(use_version_simple)
#
# use_version <- paste0(getwd(), "/src/mackerel_mvt_model_RE")
# compile(paste0(use_version, ".cpp"))
# dyn.load(use_version)
#

cross_validation <- function(x){
  set.seed(x)

  best_mod <- 1

  N_threshold = 1
  training <- Data_mackerel_use_Ireland_select %>% group_by(Release_year) %>% slice_sample(prop=P_cross) %>% ungroup()
  testing <- Data_mackerel_use_Ireland_select[-training$ID,]

  #m1 <- lm(log_rate ~ factor(Tag_area)*factor(Release_timing) + length_scaled + julian_recapture_scaled, data=training)
  m1 <- lm(log_rate ~ factor(Release_timing) + la_scaled:factor(Release_timing) + la2_scaled:factor(Release_timing) + length_scaled + julian_recapture_scaled, data=training)
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
                   N_year = length(unique(training$Release_year)),
                   Year_ID = as.numeric(as.character(training$Release_year))-2014,
                   Year_ID_pred = as.numeric(as.character(testing$Release_year))-2014,
                   Likconfig = 0)      # 0 = dnorm, 1 = dgamma
  parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                         log_sigma = rep(log(0.2),N_threshold),
                         year = rep(0, length(unique(training$Release_year))),
                         log_sigma_year = log(0.2))
  Map = list()
  Map$beta <- factor(Maps_nothres[best_mod,])
  Map$log_sigma <- factor(c(2*ncol(XX)+1))
  Map$year <- factor(ncol(XX)+1+1:length(parameters_tmb$year))
  Map$log_sigma_year <- factor(ncol(XX)+8)
  obj <- MakeADFun(data_tmb, parameters_tmb, random = "year", DLL = "mackerel_mvt_model_nothresh_RE", map=Map, silent =TRUE)
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
                    N_year = length(unique(training$Release_year)),
                    Year_ID = as.numeric(as.character(training$Release_year))-2014,
                    Year_ID_pred = as.numeric(as.character(testing$Release_year))-2014,
                    Likconfig = 0)      # 0 = dnorm, 1 = dgamma

  parameters_tmb1 <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                          log_sigma = rep(log(0.2),N_threshold),
                          year = matrix(rep(0, length(unique(training$Release_year))*N_threshold), ncol=2),
                          log_sigma_year = rep(log(0.2),N_threshold))
  Map1 = list()
  Map1$beta <- factor(Maps_best[best_mod,])
  Map1$log_sigma <- factor(c(2*ncol(XX)+1:2))
  Map1$year <- factor(N_threshold*ncol(XX)+2+matrix(1:length(parameters_tmb1$year), ncol=N_threshold))
  Map1$log_sigma_year <- factor(c(N_threshold*ncol(XX)+15:16))
  parameters_tmb1$beta[which(is.na(Maps_best[best_mod,]))] <- 0

  obj1 <- MakeADFun(data_tmb1, parameters_tmb1, random = "year", DLL = "mackerel_mvt_model_RE", map=Map1, silent =TRUE)
  opt1 <- fit_tmb( obj=obj1, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000), quiet=TRUE)
  sd_report1 <- sdreport(obj1)
  mu_pred <- matrix(summary(sd_report1, "report")[which(rownames(summary(sd_report_1break_best, "report")) == "mu_pred"),1], ncol=2, byrow=FALSE)

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

if(!do_parallel)
{
  t1  <- Sys.time()
  Cross_val_res <- 1:N_cross %>% map(function(x) cross_validation(x))
  t2  <- Sys.time()
  cat("Completed in ", t2-t1,"\n")
  Cross_val_res_df <- do.call(rbind, Cross_val_res)
  colnames(Cross_val_res_df) <- c("No_thres", "Thresh")
  boxplot(Cross_val_res_df)
  Cross_val_res_df_melt <- reshape2::melt(Cross_val_res_df, variable.name="Model_type")
  colnames(Cross_val_res_df_melt) <- c("ID", "Model_type", "MSE")
}
# ------------------------
# -- Parallell version: --
#--------------------------
if(do_parallel)
{
  library(doParallel)
  cl <- makeCluster(detectCores()-2)

  clusterEvalQ(cl, {
    library(TMB)
    library(TMBhelper)
    library(tidyverse)
    use_version_simple <- paste0(getwd(), "/src/mackerel_mvt_model_nothresh_RE")
    #compile(paste0(use_version_simple, ".cpp"))
    dyn.load(use_version_simple)
    use_version <- paste0(getwd(), "/src/mackerel_mvt_model_RE")
    #compile(paste0(use_version, ".cpp"))
    dyn.load(use_version)
  })
  clusterExport(cl,
                varlist = c("P_cross", "Data_mackerel_use_Ireland_select",
                            "threshold_vals", "mean_diff_tag_area",
                            "Maps_best", "Maps_nothres"))
  t1.core <- Sys.time()
  Cross_val_res <- t(parSapply(cl, 1:N_cross, cross_validation))
  t2.core  <- Sys.time()
  stopCluster(cl)
  cat("Completed in ", t2.core-t1.core,"\n")
  Cross_val_res_df <- data.frame(Cross_val_res)
  colnames(Cross_val_res_df) <- c("No_thres", "Thresh")
  Cross_val_res_df_melt <- gather(Cross_val_res_df, "Model_type","MSE")
}

Cross_val_res_df_melt$Model_type <- factor(Cross_val_res_df_melt$Model_type, labels=c("No threshold", "With threshold"))

Cross_val_res_df_melt %>% group_by(Model_type) %>% summarize(mean(MSE))
#expression1 <- bquote(delta[AIC] == ~ .(trunc(AIC_nothres_best-AIC_best)))
#expression2 <- bquote(delta[AIC]==0)


print(
  g1 <- ggplot(Cross_val_res_df_melt, aes(x=Model_type, y=MSE)) + geom_violin() + geom_boxplot(width=0.1, fill="lightgreen") +
    theme_bw() +
    geom_text(data=data.frame(Model_type=c("No threshold"), MSE=c(0.031)), label=deparse(expression1), parse=TRUE) +
    geom_text(data=data.frame(Model_type=c("With threshold"), MSE=c(0.031)), label=deparse(expression2), parse=TRUE) +
    coord_cartesian(ylim=c(0.017,0.032))
  g1
)

print(
  t.test(x=Cross_val_res_df[,1], y=Cross_val_res_df[,2], alternative = "two.sided", paired = TRUE)
)
# Adding a plot illustrating that the it is not just the distrbution (boxplot), but with one exception
#
print(
  g2 <- ggplot(data = as.data.frame(Cross_val_res_df), aes(x = No_thres, y = Thresh,
                                                           col = No_thres < Thresh))+
    geom_point()+
    geom_abline(slope =1, intercept = 0, lty = 2, col = 1)+
    scale_color_manual(values = c("lightblue", "red"))+
    scale_x_continuous(name = "No threshold", breaks = seq(0.0175,0.03, 0.0025), limit = range(Cross_val_res_df))+
    scale_y_continuous(name = "With threshold", breaks = seq(0.0175,0.03, 0.0025), limit = range(Cross_val_res_df))+
    theme_bw()+theme(
      panel.grid = element_blank(),
      legend.position = "none"
    )
  g2
)

grid.arrange(g1,g2,ncol=2)


