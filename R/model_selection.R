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

# option to choose the  criteria for the downward model selection
mod_sel_option <- 1     # 1. based on p-value of estimated param # 2. based on p-value on the difference of parameter estimate before ad after threshold

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

    AIC_selection[1] <- opt1break$objective * 2 + p * length(opt1break$par)
    AIC_selection2[1] <- opt1break$objective * 2 + p * length(opt1break$par)

    #-- which parameter estimate has the lowest p-value? downward model selection --
      # only removing 1
      mle <- matrix(sapply(1:length(Map$beta), function(x) ifelse(is.na(x)==F, opt1break$par[Map$beta[x]], NA)), byrow=F, ncol=2)
      cov <- solve(obj1break$he(opt1break$par))
      sde <- sqrt(diag(cov))
      pvalue_diff <- function(i){
        if (! TRUE %in% mle[i,]) out <- 2*pnorm(abs(diff(mle[i,])), 0, sqrt(cov[i,i]+cov[i+ncol(XX),i+ncol(XX)]-2*cov[i,i+ncol(XX)]), lower.tail = FALSE)
        if (TRUE %in% mle[i,]) out <- NA
        return(out)
      }
      if(mod_sel_option == 1) j <- which.max(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE)[-(length(opt1break$par)-0:1)])
      if(mod_sel_option == 2) j <- which.max(sapply(1:ncol(XX), pvalue_diff))
      Maps[(k+1):nrow(Maps), j] <- NA_real_
      parameters_tmb$beta[ifelse(j<=ncol(XX),j,j-ncol(XX)),ifelse(j/ncol(XX)>1,2,1)] <- 0
      # now fixing before and after threshold
      if (j[1] <= ncol(XX)) {Maps2[(k+1), c(j,j+ncol(XX))] <- Maps2[k,j] } #; Maps2[(k+1), (j+1):ncol(Maps2)] <- Map2$beta[(j+1):ncol(Maps2)]-1 }
      if (j[1] >  ncol(XX)) {Maps2[(k+1), c(j-ncol(XX),j)] <- Maps2[k,j-ncol(XX)]   } #; Maps2[(k+1), (j+2):ncol(Maps2)] <- Map2$beta[(j+2):ncol(Maps2)]-1 }

    # -- Now we choose the best out of these two using AIC
      Map = list()
      Map$beta <- factor(Maps[k+1,])
      Map$log_sigma <- factor(c(2*ncol(XX)+1:2))
      obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
      opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                              control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
      AIC_selection[k+1] <- opt1break$objective * 2 + p * length(opt1break$par)
      Map2 = list()
      Map2$beta <- factor(Maps2[k+1,])
      Map2$log_sigma <- factor(c(2*ncol(XX)+1:2))
      obj1break2 <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map2)
      opt1break2 <- fit_tmb( obj=obj1break2, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                              control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
      AIC_selection2[k+1] <- opt1break2$objective * 2 + p * length(opt1break2$par)

      if (AIC_selection2[k+1] < AIC_selection[k+1]) {
        Maps_best[k+1,] = Maps2[k+1,]
        parameters_tmb_best <- parameters_tmb2
      }
      if (AIC_selection[k+1] < AIC_selection2[k+1]){
        Maps_best[k+1,] = Maps[k+1,]
        parameters_tmb_best <- parameters_tmb
      }

    }
  if (k>1){ # from this step,we evaluate whether it is better to simply remove the parameter or have a single parameter for before and after threhold
    # run the best model (it is a repetition of the step above but oh well)
      Map = list()
      Map$beta <- factor(Maps_best[k,])
      Map$log_sigma <- factor(c(2*ncol(XX)+1:2))
      parameters_tmb <- parameters_tmb_best

      #-- optimize-- :
      obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
      opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                            control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

      map <- as.numeric(c(as.character(Map$beta), (length(opt1break$par)-c(1,0))))
      # which parameter estimate has the lowest p-value?
      mle <- matrix(sapply(1:length(Map$beta), function(x) ifelse(is.na(x)==F, opt1break$par[Map$beta[x]], NA)), byrow=F, ncol=2)
      cov <- solve(obj1break$he(opt1break$par))
      sde <- sqrt(diag(cov))

      if(mod_sel_option == 1) j <- which.max(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE)[-(length(opt1break$par)-0:1)])
      if(mod_sel_option == 2) j <- which.max(sapply(1:ncol(XX), pvalue_diff))
      # Correct for already removed parameters:
      j.corrected <- as.numeric(map[!is.na(map) & !duplicated(map)][j])

      # testing to remove the single parameter
        Maps[(k+1), ] <- Maps_best[k,]
        Maps[(k+1), j.corrected] <- NA_real_
        parameters_tmb$beta[ifelse(j.corrected<=ncol(XX),j.corrected,j.corrected-ncol(XX)),ifelse(j.corrected/ncol(XX)>1,2,1)] <- 0

      # testing to combine the parameter before and after the threshold
        Maps2[(k+1), ] <- Maps_best[k,]
        if (j.corrected[1] <= ncol(XX)) {Maps2[(k+1), c(j.corrected,j.corrected+ncol(XX))] <- Maps_best[k,j.corrected] } #; Maps2[(k+1), (j.corrected+1):ncol(Maps2)] <- Map2$beta[(j.corrected+1):ncol(Maps2)]-1 }
        if (j.corrected[1] >  ncol(XX)) {Maps2[(k+1), c(j.corrected-ncol(XX),j.corrected)] <- Maps_best[k,j.corrected-ncol(XX)]   } #; Maps2[(k+1), (j.corrected+2):ncol(Maps2)] <- Map2$beta[(j.corrected+2):ncol(Maps2)]-1 }

      # And now the model selection
        Map = list()
        Map$beta <- factor(Maps[k+1,])
        Map$log_sigma <- factor(c(2*ncol(XX)+1:2))
        obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
        opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                              control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
        AIC_selection[k+1] <- opt1break$objective * 2 + p * length(opt1break$par)
        Map2 = list()
        Map2$beta <- factor(Maps2[k+1,])
        Map2$log_sigma <- factor(c(2*ncol(XX)+1:2))
        obj1break2 <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map2)
        opt1break2 <- fit_tmb( obj=obj1break2, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                               control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
        AIC_selection2[k+1] <- opt1break2$objective * 2 + p * length(opt1break2$par)

        if (AIC_selection2[k+1] < AIC_selection[k+1]) {
          Maps_best[k+1,] = Maps2[k+1,]
          parameters_tmb_best <- parameters_tmb2
        }
        if (AIC_selection[k+1] <= AIC_selection2[k+1]){
          Maps_best[k+1,] = Maps[k+1,]
          parameters_tmb_best <- parameters_tmb
        }
  }

}


## Extract the best model
AIC_selection
AIC_selection2
best_mod <- c(which.min(AIC_selection),which.min(AIC_selection2))[which.min(c(AIC_selection[which.min(AIC_selection)],AIC_selection2[which.min(AIC_selection2)]))]

Map = list()
Map$beta <- factor(Maps_best[best_mod,])
Map$log_sigma <- factor(c(2*ncol(XX)+1:2))
parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                       log_sigma = rep(log(0.2),N_threshold))
parameters_tmb$beta[which(is.na(Maps_best[best_mod,]))] <- 0

#-- optimize-- :
obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                      control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))
opt1break$objective * 2 + p * length(opt1break$par)

map <- as.factor(as.numeric(c(as.character(Map$beta))))
map <- as.numeric(factor(map, labels=1:length(levels(map))))
mle <- matrix(sapply(1:length(map), function(x) ifelse(is.na(x)==F, opt1break$par[map[x]], NA)), byrow=F, ncol=2)
sde <- matrix(sapply(1:length(map), function(x) ifelse(is.na(x)==F, sqrt(diag(solve(obj1break$he(opt1break$par))))[map[x]], NA)), byrow=F, ncol=2)


cov <- solve(obj1break$he(opt1break$par))
par(mfrow=c(4,3))
for (i in 1:ncol(XX)){
if (! TRUE %in% is.na(mle[i,])) plot(density(rnorm(100000, diff(mle[i,]), sqrt(cov[i,i]+cov[i+ncol(XX),i+ncol(XX)]-2*cov[i,i+ncol(XX)])))); abline(v=0)
}

# -- alternative: --
pp <- numeric(ncol(XX))
for(i in 1:ncol(XX)){
  pp[i] <- 2*pnorm(abs(mle[i,1]-mle[i,2]), mean = 0, sd = sqrt(cov[i,i]+cov[ncol(XX)+i,ncol(XX)+i]-2*cov[i,ncol(XX)+i]), lower.tail = FALSE)
}
pp

res <- data.frame(MLE=mle[1:ncol(XX),], p_value=pp)
colnames(res) <- c("mu_before", "mu_after", "p_value")
res

plot(24-0:8, AIC_selection[1:9], xlab = "Number of beta parameters", ylab = "AIC")
abline(h=min(AIC_selection), col = 2, lty=2)
opt1break$par
which.min(AIC_selection)


