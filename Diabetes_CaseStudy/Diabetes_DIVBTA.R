library(robumeta);  library(dplyr)
options(dplyr.summarise.inform=F)
lnCVR <- function(dat,r_CT=0.6){
  # note if dependent ( abs(crossoverrho1) > 0 and < 1); note in this case dat$nT = dat$nC
  dat$lnCVR = ifelse(dat$crossoverrho1 %in% c(0,1),
                     log((dat$sT/dat$xT)/(dat$sC/dat$xC))+
                       (1/2)*(1/(dat$nT-1)-1/(dat$nC-1))+
                       (1/2)*((dat$sC^2 /(dat$nC*dat$xC^2))-(dat$sT^2 /(dat$nT*dat$xT^2))),
                     log((dat$sT/dat$xT)/(dat$sC/dat$xC))+
                       (1/2)*((dat$sC^2 /(dat$nC*dat$xC^2))-(dat$sT^2 /(dat$nT*dat$xT^2))))
  
  dat$var_lnCVR =ifelse(dat$crossoverrho1 %in% c(0,1),
                        dat$sC^2/(dat$nC*dat$xC^2) + 
                          dat$sC^4/(2*dat$nC^2*dat$xC^4)+
                          0.5*dat$nC/((dat$nC-1)^2)+0.5*dat$nT/((dat$nT-1)^2)+
                          dat$sT^2/(dat$nT*dat$xT^2) + 
                          dat$sT^4/(2*dat$nT^2*dat$xT^4),
                        dat$sC^2/(dat$nC*dat$xC^2) + dat$sT^2/(dat$nT*dat$xT^2) -
                          r_CT * 2*dat$sC*dat$sT / (dat$nC*dat$xC*dat$xT) +
                          1/(dat$nC-1) - r_CT^2/(dat$nC-1) )

  dat$lnCVR = ifelse(dat$post_SD1 %in% c(1,2),dat$lnCVR,NA )
  dat$var_lnCVR = ifelse(dat$post_SD1 %in% c(1,2),dat$var_lnCVR,NA )
  
    dat$sd_lnCVR = sqrt(dat$var_lnCVR)
    dat$p_lnCVR_normality = 2*exp(pnorm(abs(dat$lnCVR/dat$sd_lnCVR),lower.tail = F,log.p = T))
    dat$p_lnCVR_ttest = ifelse(dat$crossoverrho1 %in% c(0,1),
                               2*exp(pt(abs(dat$lnCVR/dat$sd_lnCVR),df=dat$nT+dat$nC-2,lower.tail = F,log.p = T)),
                               2*exp(pt(abs(dat$lnCVR/dat$sd_lnCVR),df=dat$nT-1,lower.tail = F,log.p = T)))
    return(dat)
}


hba1c = read.csv("github_LnCVR.csv")
hba1c$reliable = ifelse(hba1c$P_onesided < 0.025 | hba1c$P_onesided > 0.975 | is.na(hba1c$P_onesided),
                        0,1)
hba1c = lnCVR(hba1c) # computes lnCVR (also included in file)

# 3 analyses
# 1) Use all 227 unique trials 
# 2) restrict to 152 reliable trials when deciding outlier criteria; reliable==1
# Definition of reliable 
#   Restrict to those for which we can compute an lnCVR
#   Exclude trials with an improbable baseline data distribution (one-sided p < 0.025 or > 0.975)
#   Exclude small trials (post-trial n1+n2 < 30)
#   Exclude cross-over trials
# 3) Bootstrapping 

# Non-parametric; use Tukey fences
nonparametric_approach <- function(dataset,reliable_values){
  # take the median across lnCVR for trials with multiple values
  hba1c_median_acrossPMID =  dataset %>% filter(!is.na(lnCVR)) %>% group_by(PMID) %>% 
    summarize(lnCVRmedian = median(lnCVR),sig = ifelse(sum(p_lnCVR_ttest < 0.05),1,0))
  # using the medians -> compute quartiles + IQR to define outliers
  hba1c_stat = dataset %>% filter(!is.na(lnCVR) & reliable %in% reliable_values) %>% group_by(PMID) %>% 
    summarize(lnCVRmedian = median(lnCVR)) %>%  ungroup() %>%  summarize(Q1 = quantile(lnCVRmedian, 0.25,na.rm = T),
                             Q3 = quantile(lnCVRmedian, 0.75,na.rm = T), 
                             IQR = IQR(lnCVRmedian,na.rm = T)) %>%
    mutate(outlier_L = Q1 - 1.5 * IQR, outlier_U = Q3 + 1.5 * IQR,
           faroutlier_L = Q1 - 3 * IQR, faroutlier_U = Q3 + 3 * IQR)
  
  hba1c_median_acrossPMID$TukeyOutlier = ifelse(hba1c_median_acrossPMID$lnCVRmedian > hba1c_stat$faroutlier_U |
                                                  hba1c_median_acrossPMID$lnCVRmedian < hba1c_stat$faroutlier_L,
                                                "Far outlier",
                                                ifelse(hba1c_median_acrossPMID$lnCVRmedian > hba1c_stat$outlier_U | 
                                                         hba1c_median_acrossPMID$lnCVRmedian < hba1c_stat$outlier_L,
                                                       "Outlier","Not outlier"))
  return(list(hba1c_stat,hba1c_median_acrossPMID))
}

# Parametric approach; use all trials and use robumeta to combine across lnCVR
# Helper function to calculate prediction intervals
robu_predict_interval <- function(fit, level = 0.95) {
  # For df we can use N-2 as suggested by Introduction to Meta-Analysis (2nd edition)
  # Chapter 17 page 122 or fit$dfs provided by robumeta
  df_val = fit$N - 2 ; p_crit <- (1 - level) / 2
  return(as.numeric(fit$reg_table[2]) + c(-1, 1) *
           qt(p = p_crit, df = df_val, lower.tail = FALSE) *
           sqrt(as.numeric(fit$mod_info$tau.sq) + as.numeric(fit$reg_table[3])^2))
}

parametric_approach <- function(dataset,reliable_values){
  hba1c_nonNA =  dataset %>% filter(!is.na(lnCVR)) 
  metaanalysis_lnCVRs <- robu(formula = lnCVR ~  1, data = hba1c_nonNA[hba1c_nonNA$reliable %in% reliable_values,],
                              studynum = PMID, var.eff.size = var_lnCVR, rho = 0.8,small = TRUE)       
  sigma_lvl_diabetes = as.data.frame(cbind(Sigma = c(3,4,5,6),
                                           rbind(robu_predict_interval(metaanalysis_lnCVRs, level = 1-(1-pnorm(3))*2),
                                                 robu_predict_interval(metaanalysis_lnCVRs, level = 1-(1-pnorm(4))*2),
                                                 robu_predict_interval(metaanalysis_lnCVRs, level = 1-(1-pnorm(5))*2),
                                                 robu_predict_interval(metaanalysis_lnCVRs, level = 1-(1-pnorm(6))*2))))
  colnames(sigma_lvl_diabetes) = c("Sigma","L","U")

  hba1c_nonNA$sig = ifelse(hba1c_nonNA$p_lnCVR_ttest < 0.05 ,1,0)
  hba1c_nonNA$Outlier = ifelse(hba1c_nonNA$lnCVR < sigma_lvl_diabetes$L[sigma_lvl_diabetes$Sigma == 6] |
                                 hba1c_nonNA$lnCVR > sigma_lvl_diabetes$U[sigma_lvl_diabetes$Sigma == 6],
                               "6+",
                               ifelse(hba1c_nonNA$lnCVR < sigma_lvl_diabetes$L[sigma_lvl_diabetes$Sigma == 5] |
                                        hba1c_nonNA$lnCVR > sigma_lvl_diabetes$U[sigma_lvl_diabetes$Sigma == 5],
                                      "5",
                                      ifelse(hba1c_nonNA$lnCVR < sigma_lvl_diabetes$L[sigma_lvl_diabetes$Sigma == 4] |
                                               hba1c_nonNA$lnCVR > sigma_lvl_diabetes$U[sigma_lvl_diabetes$Sigma == 4],
                                             "4",
                                             ifelse(hba1c_nonNA$lnCVR < sigma_lvl_diabetes$L[sigma_lvl_diabetes$Sigma == 3] |
                                                      hba1c_nonNA$lnCVR > sigma_lvl_diabetes$U[sigma_lvl_diabetes$Sigma == 3],
                                                    "3","Not outlier"))))
  
  return(list(sigma_lvl_diabetes,hba1c_nonNA))
}


bootstrap_helperfcn <- function(dataset,reliable_values=c(0,1)){
  pmid_sample = sample(unique(dataset$PMID),size=length(unique(dataset$PMID)),replace = T)
  bootstrap_sample=NULL
  for(i in 1:length(pmid_sample)){
    add = dataset[dataset$PMID == pmid_sample[i],]; add$newStudyNum = i
    bootstrap_sample = rbind(bootstrap_sample,add)
  }
  bootstrap_sample$PMID = bootstrap_sample$newStudyNum
  nonparametric_boot = nonparametric_approach(dataset = bootstrap_sample,reliable_values=c(0,1))
  parametric_boot = parametric_approach(dataset = bootstrap_sample,reliable_values=c(0,1))
  return(list(nonparametric_boot[[1]],parametric_boot[[1]]))
}
  
bootstrap_fcn <- function(dataset,reliable_values=c(0,1),initial_seed=1127,bootstrap_samples){
  set.seed(initial_seed); nonparametric_boot_df = parametric_boot_df = NULL
  for(j in 1:bootstrap_samples){
    boot_samp = bootstrap_helperfcn(dataset=dataset,reliable_values=reliable_values)
    nonparametric_boot_samp = as.data.frame(boot_samp[[1]]); nonparametric_boot_samp$Boot = j
    parametric_boot_samp = as.data.frame(boot_samp[[2]]); parametric_boot_samp$Boot = j
    nonparametric_boot_df = rbind(nonparametric_boot_df,nonparametric_boot_samp)
    parametric_boot_df = rbind(parametric_boot_df,parametric_boot_samp)
  }
  # compute 95% CI for all variables of interest
  nonparametric_boot_vals = nonparametric_boot_df %>% reframe(
    value = c("far outlier L","outlier L","outlier U","far outlier U"),
    L = c(quantile(faroutlier_L,0.025),quantile(outlier_L,0.025),
          quantile(outlier_U,0.025),quantile(faroutlier_U,0.025)),
    Est = c(quantile(faroutlier_L,0.5),quantile(outlier_L,0.5),
          quantile(outlier_U,0.5),quantile(faroutlier_U,0.5)),
    U = c(quantile(faroutlier_L,0.975),quantile(outlier_L,0.975),
          quantile(outlier_U,0.975),quantile(faroutlier_U,0.975)))
  
  parametric_boot_vals = parametric_boot_df %>% group_by(Sigma) %>% 
    summarize(L_boot_L =  quantile(L,0.025),L_boot = quantile(L,0.5),L_boot_U = quantile(L,0.975),
              U_boot_L =  quantile(U,0.025),U_boot = quantile(U,0.5),U_boot_U = quantile(U,0.975))
    
  return(list(nonparametric_boot_vals,parametric_boot_vals))
}
# info:
length(unique(hba1c$PMID[!is.na(hba1c$lnCVR)])); sum(!is.na(hba1c$lnCVR))
length(unique(hba1c$PMID[!is.na(hba1c$lnCVR) & hba1c$reliable==1])); sum(!is.na(hba1c$lnCVR) & hba1c$reliable==1)

  
# 1) Use all 227 unique trials 
nonparametric_all = nonparametric_approach(dataset = hba1c,reliable_values=c(0,1))
nonparametric_all_df = nonparametric_all[[2]]
parametric_all = parametric_approach(dataset = hba1c,reliable_values=c(0,1))
parametric_all_df = parametric_all[[2]]; 
# 2) restrict to 152 reliable trials when deciding outlier criteria
nonparametric_reliable = nonparametric_approach(hba1c,reliable_values=c(1))
nonparametric_reliable_df = nonparametric_reliable[[2]]
parametric_reliable = parametric_approach(hba1c,reliable_values=c(1))
parametric_reliable_df = parametric_reliable[[2]]; 
# prediction intervals:
parametric_reliable[[1]]; parametric_all[[1]]
nonparametric_reliable[[1]]; nonparametric_all[[1]]
# identify significant lnCVRs:
nrow(nonparametric_all_df[nonparametric_all_df$sig ==1,])
length(unique(nonparametric_all_df$PMID[nonparametric_all_df$sig ==1]))

nrow(parametric_all_df[parametric_all_df$sig ==1,])
length(unique(parametric_all_df$PMID[parametric_all_df$sig ==1]))

# identify significant outliers:
nrow(nonparametric_all_df[nonparametric_all_df$sig ==1 & nonparametric_all_df$TukeyOutlier != "Not outlier",])
length(unique(nonparametric_all_df$PMID[nonparametric_all_df$sig ==1 & nonparametric_all_df$TukeyOutlier != "Not outlier"]))

nrow(nonparametric_reliable_df[nonparametric_reliable_df$sig ==1 & nonparametric_reliable_df$TukeyOutlier != "Not outlier",])
length(unique(nonparametric_reliable_df$PMID[nonparametric_reliable_df$sig ==1 & nonparametric_reliable_df$TukeyOutlier != "Not outlier"]))

nrow(parametric_all_df[parametric_all_df$sig ==1 & parametric_all_df$Outlier != "Not outlier",])
length(unique(parametric_all_df$PMID[parametric_all_df$sig ==1 & parametric_all_df$Outlier != "Not outlier"]))

nrow(parametric_reliable_df[parametric_reliable_df$sig ==1 & parametric_reliable_df$Outlier != "Not outlier",])
length(unique(parametric_reliable_df$PMID[parametric_reliable_df$sig ==1 & parametric_reliable_df$Outlier != "Not outlier"]))

# 3) Bootstrapping 
bootstrap_vals = bootstrap_fcn(dataset=hba1c,reliable_values = c(0,1),bootstrap_samples = 1000)
bootstrap_vals

# 4) Do any of of the identified problematic trials have significant outlier effect sizes 
metaanalysis_effectsizes_lnCVRs <- robu(formula = effect ~  1, data = hba1c[!is.na(hba1c$effect),],
                            studynum = PMID, var.eff.size = var_effect, rho = 0.8,small = TRUE)       
sigma_lvl_diabetes_effectsizes = as.data.frame(cbind(Sigma = c(3,4,5,6),
                                         rbind(robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(3))*2),
                                               robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(4))*2),
                                               robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(5))*2),
                                               robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(6))*2))))
colnames(sigma_lvl_diabetes_effectsizes) = c("Sigma","L","U")
sigma_lvl_diabetes_effectsizes

# any of detected identified problematic trials have significant outlier effect sizes :
parametric_reliable_df[parametric_reliable_df$sig ==1 & parametric_reliable_df$Outlier != "Not outlier" &
                        abs(parametric_reliable_df$effect/sqrt(parametric_reliable_df$var_effect)) > abs(qnorm(p=0.025)) & 
                         (parametric_reliable_df$effect <  sigma_lvl_diabetes_effectsizes$L[sigma_lvl_diabetes_effectsizes$Sigma==3] |
                            parametric_reliable_df$effect > sigma_lvl_diabetes_effectsizes$U[sigma_lvl_diabetes_effectsizes$Sigma==3]),]
