# -----------------
# -- Full model: --
# -----------------
Map = list(beta = matrix(1:(ncol(XX)*2), ncol = 2),
           log_sigma = c(2*ncol(XX)+1:2))
Map$beta <- factor(Map$beta)
Map$log_sigma <- factor(Map$log_sigma)

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
                 Likconfig = 0      # 0 = dnorm, 1 = dgamma
)

parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                       log_sigma = rep(log(0.2),N_threshold)
)
data_tmb$Likconfig = 0
obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE, control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
opt1break
sde <- sqrt(diag(solve(obj1break$he(opt1break$par))))
sde_beta <- matrix(sde[1:(length(opt1break$par)-2)], ncol = 2)
matrix(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE), ncol = 2)

# ---------------------
# -- Model selection --
# ---------------------

# Maps selects which parameters are set to zero, each row representing one model:
Maps <- matrix(1:(ncol(XX)*2), ncol = (ncol(XX)*2), nrow =ncol(XX)+1,byrow=TRUE)
Maps2 <- matrix(1:(ncol(XX)*2), ncol = (ncol(XX)*2), nrow =ncol(XX)+1,byrow=TRUE)
Maps_best <- matrix(1:(ncol(XX)*2), ncol = (ncol(XX)*2), nrow =ncol(XX)+1,byrow=TRUE)
AIC_selection <- numeric(nrow(Maps)-1)
AIC_selection2 <- numeric(nrow(Maps)-1)
p = 2 # AIC
#p = log(length(data_tmb$y)) # BIC

# -- Initial values: --
parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                       log_sigma = rep(log(0.2),N_threshold))
parameters_tmb2 <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                       log_sigma = rep(log(0.2),N_threshold))
parameters_tmb_best <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                       log_sigma = rep(log(0.2),N_threshold))



# -- looping over potenial models: --
for(k in 1:(nrow(Maps)-1)){
  if (k==1){
    Map = list()
    Map$beta <- factor(Maps[k,])
    Map$log_sigma <- factor(c(2*ncol(XX)+1:2))

    #-- optimize-- :
    obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
    opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                          control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

    # -- AIC: --

    AIC_selection[k] <- opt1break$objective * 2 + p * length(opt1break$par)
    # which parameter estimate has the lowest p-value?
    sde <- sqrt(diag(solve(obj1break$he(opt1break$par))))
    j <- which.max(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE)[-(length(opt1break$par)-0:1)])
    Maps[(k+1):nrow(Maps), j] <- NA_real_
    parameters_tmb$beta[j%%ncol(XX),ifelse(j/ncol(XX)>1,2,1)] <- 0

    if (j%%2==0) {Maps2[(k+1):nrow(Maps2), (j-1):j] <- Maps2[k,j-1] ; Maps2[(k+1), (j+1):ncol(Maps2)] <- Maps2[k, (j+1):ncol(Maps2)]-1 }
    if (j%%2==1) {Maps2[(k+1):nrow(Maps2), j:(j+1)] <- Maps2[k,j]   ; Maps2[(k+1), (j+2):ncol(Maps2)] <- Maps2[k, (j+2):ncol(Maps2)]-1 }


    }
  if (k==2){ # from this step,we evaluate whether it is better to simply remove the parameter or have a single parameter for before and after threhold
    # testing to remove the single parameter
      Map = list()
      Map$beta <- factor(Maps[k,])
      Map$log_sigma <- factor(c(2*ncol(XX)+1:2))

      #-- optimize-- :
      obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
      opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                            control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

      # -- AIC: --
      map <- as.numeric(c(as.character(Map$beta), (length(opt1break$par)-c(1,0))))
      mle <- sapply(1:length(map), function(x) ifelse(is.na(x)==F, opt1break$par[map[x]], NA))
      names(mle) <- 1:length(map)
      AIC_selection[k] <- opt1break$objective * 2 + p * length(opt1break$par)
      # which parameter estimate has the lowest p-value?
      sde <- sqrt(diag(solve(obj1break$he(opt1break$par))))
      j <- which.max(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE)[-(length(opt1break$par)-0:1)])
      # Correct for already removed parameters:
      index.name <- as.numeric(names(mle[!is.na(mle)][j]))
      j.corrected <- which(Maps[k,] == index.name)
      Maps[(k+1):nrow(Maps), j.corrected] <- NA_real_
      parameters_tmb$beta[j.corrected%%ncol(XX),ifelse(j.corrected/ncol(XX)>1,2,1)] <- 0

    # testing to combine the parameter before and after the threshold
      Map2 = list()
      Map2$beta <- factor(Maps2[k,])
      Map2$log_sigma <- factor(c(2*ncol(XX)+1:2))

      #-- optimize-- :
      obj1break2 <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map2)
      opt1break2 <- fit_tmb( obj=obj1break2, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                            control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

      # -- AIC: --
      map2 <- as.numeric(c(as.character(Map2$beta), (length(opt1break2$par)-c(1,0))))
      mle2 <- sapply(1:length(map2), function(x) ifelse(is.na(x)==F, opt1break2$par[map2[x]], NA))
      names(mle2) <- 1:length(map2)

      #names(mle2) <- as.numeric(c(as.character(Map2$beta[!is.na(Map2$beta)]), as.character(Map2$log_sigma[!is.na(Map2$log_sigma)])))
      AIC_selection2[k] <- opt1break2$objective * 2 + p * length(opt1break2$par)
      # which parameter estimate has the lowest p-value?
      sde <- sqrt(diag(solve(obj1break2$he(opt1break2$par))))
      j <- which.max(2*pnorm(abs(opt1break2$par)/sde, lower.tail = FALSE)[-(length(opt1break2$par)-0:1)])
      # Correct for already removed parameters:
      index.name <- as.numeric(names(mle2[!is.na(mle2)][j]))
      j.corrected <- which(Maps2[k,] == j)

      if (j.corrected%%2==0) {Maps2[(k+1):nrow(Maps2), (j.corrected-1):j.corrected] <- Maps2[k,j.corrected-1] ; Maps2[(k+1), (j.corrected+1):ncol(Maps2)] <- Maps2[k, (j.corrected+1):ncol(Maps2)]-1 }
      if (j.corrected%%2==1) {Maps2[(k+1):nrow(Maps2), j.corrected:(j.corrected+1)] <- Maps2[k,j.corrected]   ; Maps2[(k+1), (j.corrected+2):ncol(Maps2)] <- Maps2[k, (j.corrected+2):ncol(Maps2)]-1 }

    # And now, we choose the best out of these two
      if (AIC_selection2[k] < AIC_selection[k]) {
        Maps_best[k+1,] = Maps2[k+1,]
        parameters_tmb_best <- parameters_tmb2
      }
      if (AIC_selection[k] < AIC_selection2[k]){
        Maps_best[k+1,] = Maps[k+1,]
        parameters_tmb_best <- parameters_tmb
      }
  }
    if (k>2){ # from this step,we evaluate whether it is better to simply remove the parameter or have a single parameter for before and after threhold
      # testing to remove the single parameter
      Map = list()
      Map$beta <- factor(Maps_best[k,])
      Map$log_sigma <- factor(c(2*ncol(XX)+1:2))
      parameters_tmb <- parameters_tmb_best

      #-- optimize-- :
      obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
      opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                            control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

      # -- AIC: --
      map <- as.numeric(c(as.character(Map$beta), (length(opt1break$par)-c(1,0))))
      mle <- sapply(1:length(map), function(x) ifelse(is.na(x)==F, opt1break$par[map[x]], NA))
      names(mle) <- 1:length(map)
      AIC_selection[k] <- opt1break$objective * 2 + p * length(mle)
      # which parameter estimate has the lowest p-value?
      sde <- sqrt(diag(solve(obj1break$he(opt1break$par))))
      j <- which.max(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE)[-(length(opt1break$par)-0:1)])
      # Correct for already removed parameters:
      index.name <- as.numeric(names(mle[!is.na(mle)][j]))
      j.corrected <- which(Maps[k,] == index.name)
      Maps[(k+1):nrow(Maps), j.corrected] <- NA_real_
      parameters_tmb$beta[j.corrected%%ncol(XX),ifelse(j.corrected/ncol(XX)>1,2,1)] <- 0

      # testing to combine the parameter before and after the threshold
      Map2 = list()
      Map2$beta <- factor(Maps_best[k,])
      Map2$log_sigma <- factor(c(2*ncol(XX)+1:2))
      parameters_tmb2 <- parameters_tmb_best

      #-- optimize-- :
      obj1break2 <- MakeADFun(data_tmb, parameters_tmb2, random = NULL, DLL = "mackerel_mvt_model", map=Map2)
      opt1break2 <- fit_tmb( obj=obj1break2, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                             control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

      # -- AIC: --
      map2 <- as.numeric(c(as.character(Map2$beta), (length(opt1break2$par)-c(1,0))))
      mle2 <- sapply(1:length(map2), function(x) ifelse(is.na(x)==F, opt1break2$par[map2[x]], NA))
      names(mle2) <- 1:length(map2)

      #names(mle2) <- as.numeric(c(as.character(Map2$beta[!is.na(Map2$beta)]), as.character(Map2$log_sigma[!is.na(Map2$log_sigma)])))
      AIC_selection2[k] <- opt1break2$objective * 2 + p * length(opt1break2$par)
      # which parameter estimate has the lowest p-value?
      sde <- sqrt(diag(solve(obj1break2$he(opt1break2$par))))
      j <- which.max(2*pnorm(abs(opt1break2$par)/sde, lower.tail = FALSE)[-(length(opt1break2$par)-0:1)])
      # Correct for already removed parameters:
      index.name <- as.numeric(names(mle2[!is.na(mle2)][j]))
      j.corrected <- which(Map2$beta == j)

      if (j.corrected%%2==0) {Maps2[(k+1):nrow(Maps2), (j.corrected-1):j.corrected] <- Maps2[k,j.corrected-1] ; Maps2[(k+1), (j.corrected+1):ncol(Maps2)] <- Maps2[k, (j.corrected+1):ncol(Maps2)]-1 }
      if (j.corrected%%2==1) {Maps2[(k+1):nrow(Maps2), j.corrected:(j.corrected+1)] <- Maps2[k,j.corrected]   ; Maps2[(k+1), (j.corrected+2):ncol(Maps2)] <- Maps2[k, (j.corrected+2):ncol(Maps2)]-1 }

      # And now, we choose the best out of these two
      if (AIC_selection2[k] < AIC_selection[k]) {
        Maps_best[k+1,] = Maps2[k+1,]
        parameters_tmb_best <- parameters_tmb2
      }
      if (AIC_selection[k] < AIC_selection2[k]){
        Maps_best[k+1,] = Maps[k+1,]
        parameters_tmb_best <- parameters_tmb
      }

    }

}



plot(24-0:8, AIC_selection[1:9], xlab = "Number of beta parameters", ylab = "AIC")
abline(h=min(AIC_selection), col = 2, lty=2)
opt1break$par
which.min(AIC_selection)
