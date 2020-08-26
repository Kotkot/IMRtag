# -----------------
# -- Full model: --
# -----------------
Map = list(beta = matrix(1:22, ncol = 2),
           log_sigma = c(23:24))
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
Maps <- matrix(1:22, ncol = 22, nrow =12+1,byrow=TRUE)
AIC_selection <- numeric(nrow(Maps)-1)
p = 2 # AIC
#p = log(length(data_tmb$y)) # BIC

# -- Initial values: --
parameters_tmb <- list(beta = matrix(c(rep(10,N_threshold),runif((ncol(XX)-1)*N_threshold,-2,2)),byrow=T, ncol=N_threshold),
                       log_sigma = rep(log(0.2),N_threshold)
)


# -- looping over potenial models: --
for(k in 1:(nrow(Maps)-1)){
Map = list()
Map$beta <- factor(Maps[k,])
Map$log_sigma <- factor(25:26)

#-- optimize-- :
obj1break <- MakeADFun(data_tmb, parameters_tmb, random = NULL, DLL = "mackerel_mvt_model", map=Map)
opt1break <- fit_tmb( obj=obj1break, lower=-14, upper=14, getsd=FALSE, bias.correct=FALSE,
                      control = list(eval.max = 20000, iter.max = 20000, trace = TRUE))

# -- AIC: --

AIC_selection[k] <- opt1break$objective * 2 + p * length(opt1break$par)
# which parameter estimate has the lowest p-value?
sde <- sqrt(diag(solve(obj1break$he(opt1break$par))))
j <- which.max(2*pnorm(abs(opt1break$par)/sde, lower.tail = FALSE)[-(length(opt1break$par)-0:1)])
# Correct for already removed parameters:
j.corrected <- (Maps[k,!is.na(Maps[k,])])[j]
Maps[(k+1):nrow(Maps), j.corrected] <- NA_real_
parameters_tmb$beta[j.corrected] <- 0
}
plot(24-0:8, AIC_selection[1:9], xlab = "Number of beta parameters", ylab = "AIC")
abline(h=min(AIC_selection), col = 2, lty=2)
opt1break$par
which.min(AIC_selection)
