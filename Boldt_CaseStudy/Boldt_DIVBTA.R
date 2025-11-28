library(robumeta);  library(dplyr)

diff_var <- function(dat){
  # ratio of variances 
  dat$ratio_var = dat$sT^2/dat$sC^2
  dat$p_ratio_var = 2 * ifelse(dat$ratio_var < 1,
                           pf(dat$ratio_var, df1=dat$nT-1, df2=dat$nC-1, lower.tail = TRUE),
                           pf(dat$ratio_var, df1=dat$nT-1, df2=dat$nC-1, lower.tail = FALSE))
  # log of the ratio of SDs (lnVR)
  dat$lnVR = log(dat$sT/dat$sC)+1/(2*(dat$nT-1))-1/(2*(dat$nC-1))
  dat$sd_lnVR = sqrt(1/(2*(dat$nT-1))+1/(2*(dat$nC-1)))
  dat$var_lnVR =  dat$sd_lnVR^2
  dat$p_lnVR_norm = 2*exp(pnorm(abs(dat$lnVR/dat$sd_lnVR),lower.tail = F,log.p = T))
  dat$p_lnVR_t = 2*exp(pt(abs(dat$lnVR/dat$sd_lnVR),df = dat$nT + dat$nC - 2,lower.tail = F,log.p = T))
  return(dat)
}

# Data from Boldt abstract:
# -5.9 +/- 1.2 mmol/L vs +0.2 +/- 0.2 mmol/L (n=25)
diff_var(data.frame(author="Boldt",
                    xT=-5.9,sT=1.2,nT=25,
                    xC=0.2,sC=0.2,nC=25))[,c("p_ratio_var","p_lnVR_norm","p_lnVR_t")]

comparatorValues = read.csv("base.csv")
comparatorValues = diff_var(comparatorValues)

# 
diff_var(data.frame(author="Boldt",
                    xT=-5.9,sT=1.2,nT=25,
                    xC=0.2,sC=0.2,nC=25))[,c("lnVR")]

# Non-parametric distribution; Tukey fences:
# outlier if beyond: (Q1 - 1.5 * IQR or Q3 + 1.5 * IQR), 
# far outlier if beyond: (Q1 - 3 * IQR or Q3 + 3 * IQR).

comparatorValues %>%
  summarise(Q1 = quantile(lnVR, 0.25),Q3 = quantile(lnVR, 0.75), IQR = IQR(lnVR)) %>%
  mutate(mild_L = Q1 - 1.5 * IQR, mild_U = Q3 + 1.5 * IQR,
         extreme_L = Q1 - 3 * IQR, extreme_U = Q3 + 3 * IQR) %>%
  select(extreme_L,mild_L,mild_U,extreme_U)

# Parametric distributions 
metaanalysis_lnVRs <- robu(formula = lnVR ~  1, data = comparatorValues,
                                studynum = studynum, var.eff.size = var_lnVR, rho = 0.8,small = TRUE)       

# Function to calculate prediction intervals
robu_predict_interval <- function(fit, level = 0.95) {
  # For df we can use N-2 as suggested by Introduction to Meta-Analysis (2nd edition)
  # Chapter 17 page 122 or fit$dfs provided by robumeta
  df_val = fit$N - 2 ; p_crit <- (1 - level) / 2
  return(as.numeric(fit$reg_table[2]) + c(-1, 1) *
           qt(p = p_crit, df = df_val, lower.tail = FALSE) *
           sqrt(as.numeric(fit$mod_info$tau.sq) + as.numeric(fit$reg_table[3])^2))
}

parametric_intervals = cbind(Sigma = c(3,4,5,6),
      rbind(robu_predict_interval(metaanalysis_lnVRs, level = 1-(1-pnorm(3))*2),
            robu_predict_interval(metaanalysis_lnVRs, level = 1-(1-pnorm(4))*2),
            robu_predict_interval(metaanalysis_lnVRs, level = 1-(1-pnorm(5))*2),
            robu_predict_interval(metaanalysis_lnVRs, level = 1-(1-pnorm(6))*2)))
colnames(parametric_intervals) = c("Sigma level","L","U")
(parametric_intervals)


