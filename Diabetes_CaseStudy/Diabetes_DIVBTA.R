library(robumeta);  library(dplyr)
options(dplyr.summarise.inform=F)
lnCVR <- function(dat,r_CT=0.6){
  # note if dependent ( abs(crossoverrho1) > 0 and < 1); note in this case dat$nT = dat$nC
  # note if crossoverrho1==1 that means only first period was used
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

  dat$lnCVR = ifelse(dat$SD_Available ==1 ,dat$lnCVR,NA )
  dat$var_lnCVR = ifelse(dat$SD_Available ==1,dat$var_lnCVR,NA )

  dat$V_T = dat$sT^2/(dat$nT*dat$xT^2) + dat$sT^4/(2*dat$nT^2*dat$xT^4) + 0.5*dat$nT/((dat$nT-1)^2)
  dat$V_C = dat$sC^2/(dat$nC*dat$xC^2) +  dat$sC^4/(2*dat$nC^2*dat$xC^4)+0.5*dat$nC/((dat$nC-1)^2)
  
  # Welch-Satterthwaite df for parallel-arm lnCVR
  dat$df_welch =ifelse(dat$crossoverrho1 %in% c(0,1),
                       (dat$V_T + dat$V_C)^2 / ( (dat$V_T^2 / (dat$nT - 1)) + (dat$V_C^2 / (dat$nC - 1)) ),
                       pmin(dat$nT,dat$nC)-1)
  
    dat$sd_lnCVR = sqrt(dat$var_lnCVR)
    dat$p_lnCVR_normality = 2*exp(pnorm(abs(dat$lnCVR/dat$sd_lnCVR),lower.tail = F,log.p = T))
    dat$p_lnCVR_ttest = ifelse(dat$crossoverrho1 %in% c(0,1),
                               2*exp(pt(abs(dat$lnCVR/dat$sd_lnCVR),df=dat$nT+dat$nC-2,lower.tail = F,log.p = T)),
                               2*exp(pt(abs(dat$lnCVR/dat$sd_lnCVR),df=dat$nT-1,lower.tail = F,log.p = T)))
    dat$p_lnCVR_welch =  2*exp(pt(abs(dat$lnCVR/dat$sd_lnCVR),df=dat$df_welch,lower.tail = F,log.p = T))
    return(dat)
}


hba1c <- read_sas("plos.sas7bdat")
hba1c$reliable = ifelse(hba1c$P_onesided < 0.025 | hba1c$P_onesided > 0.975 | is.na(hba1c$P_onesided),
                        0,1)
hba1c = lnCVR(hba1c) # computes lnCVR (also included in file)

# 3 analyses
# 1) Use all 226 unique trials 
# 2) restrict to 175 reliable trials when deciding outlier criteria; reliable==1
# Definition of reliable 
#   Restrict to those for which we can compute an lnCVR
#   Exclude trials with an improbable baseline data distribution (one-sided p < 0.025 or > 0.975 or missing)
# 3) Bootstrapping 

# Non-parametric; use Tukey fences
nonparametric_approach <- function(dataset,reliable_values){
  # take the median across lnCVR for trials with multiple values
  hba1c_median_acrossPMID =  dataset %>% filter(!is.na(lnCVR)) %>% group_by(PMID) %>% 
    summarize(lnCVRmedian = median(lnCVR),sig = ifelse(sum(p_lnCVR_welch < 0.05),1,0)) 
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
  # Introduction to Meta-Analysis (2nd edition)  Chapter 17 page 122
  # For df we can use N-2 for independent observations 
  # or fit$dfs provided by robumeta for small = TRUE which is the recommended default option.
  
  df_val <-  fit$reg_table$dfs[1]  #  df_val = fit$N - 2 ;
  p_crit <- (1 - level) / 2
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

  hba1c_nonNA$sig = ifelse(hba1c_nonNA$p_lnCVR_welch < 0.05 ,1,0)
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

  
# 1) Use all 226 unique trials 
nonparametric_all = nonparametric_approach(dataset = hba1c,reliable_values=c(0,1))
nonparametric_all_df = nonparametric_all[[2]]
parametric_all = parametric_approach(dataset = hba1c,reliable_values=c(0,1))
parametric_all_df = parametric_all[[2]]; 
# 2) restrict to 175 reliable trials when deciding outlier criteria
nonparametric_reliable = nonparametric_approach(hba1c,reliable_values=c(1))
nonparametric_reliable_df = nonparametric_reliable[[2]]
parametric_reliable = parametric_approach(hba1c,reliable_values=c(1))
parametric_reliable_df = parametric_reliable[[2]]; 
# 
hba1c_nonNA_temp =  hba1c %>% filter(!is.na(lnCVR)) ;reliable_values=c(1)
(metaanalysis_lnCVRs <- robu(formula = lnCVR ~  1, data = hba1c_nonNA_temp[hba1c_nonNA_temp$reliable %in% reliable_values,],
                            studynum = PMID, var.eff.size = var_lnCVR, rho = 0.8,small = TRUE))   

# prediction intervals:
parametric_reliable[[1]]; parametric_all[[1]]

# Sigma          L         U
# 1     3 -0.5378088 0.4662464
# 2     4 -0.7135134 0.6419510
# 3     5 -0.8967035 0.8251410
# 4     6 -1.0896645 1.0181020
# Sigma          L         U
# 1     3 -0.5784954 0.5338991
# 2     4 -0.7707704 0.7261741
# 3     5 -0.9691333 0.9245370
# 4     6 -1.1753567 1.1307604

nonparametric_reliable[[1]]; nonparametric_all[[1]]

# # A tibble: 1 × 7
# Q1     Q3   IQR outlier_L outlier_U faroutlier_L faroutlier_U
# <dbl>  <dbl> <dbl>     <dbl>     <dbl>        <dbl>        <dbl>
#   1 -0.166 0.0784 0.244    -0.533     0.445       -0.899        0.812
# # A tibble: 1 × 7
# Q1     Q3   IQR outlier_L outlier_U faroutlier_L faroutlier_U
# <dbl>  <dbl> <dbl>     <dbl>     <dbl>        <dbl>        <dbl>
#   1 -0.168 0.0982 0.266    -0.567     0.497       -0.966        0.896


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
# 
print(as.data.frame(parametric_reliable_df[parametric_reliable_df$sig ==1 & parametric_reliable_df$Outlier != "Not outlier",
                                           c("PMID", "First_Author","crossoverrho1",  "lnCVR", "var_lnCVR","sigma_lncvr", "p_lnCVR_welch")]))

# PMID     First_Author crossoverrho1      lnCVR  var_lnCVR sigma_lncvr p_lnCVR_welch
# 9      -1504           Fang Q          0.00  1.0441257 0.01911739          6+  1.411316e-11
# 16     -1300          Huang X          0.00  0.5114571 0.04206064           3  1.556432e-02
# 36   3106656        Schade DS          0.03  0.5763776 0.04725572           3  1.897612e-02
# 47   6758461         Meschi F          0.00 -1.0306233 0.17405035           5  2.699940e-02
# 65   9314625 Agurs-Collins TD          0.00 -0.6399895 0.04677135           3  4.728877e-03
# 90  14513883          Cohen D          0.68  1.1326443 0.05970280          6+  7.219178e-04
# 105 15766369          Kiran M          0.00 -0.8381108 0.05426024           4  8.406900e-04
# 127 17379048         Zieve FJ          0.00 -0.6147961 0.03749808           3  2.446822e-03
# 130 17493509          Dans AM          0.00  1.4706996 0.05662549          6+  3.225277e-07
# 131 17561800         Homko CJ          0.00 -0.9977952 0.04393745           5  1.747480e-05
# 141 17909165        Vincent D          0.00 -0.8545878 0.15720772           4  4.824750e-02
# 143 18058596         Schiel R          0.00 -0.5424229 0.06804593           3  4.567815e-02
# 154 18715214        Hirsch IB          0.00  0.5213207 0.06050650           3  4.195709e-02
# 182 19792844    Al-Zahrani MS          0.00  0.9015284 0.06136737           5  1.037316e-03
# 202 20923485       Dur\xe1n A          0.00 -0.5639796 0.01378154           3  4.270340e-06
# 224 21751889      Petrovski G          0.00  0.6665609 0.09885477           4  4.515334e-02
# 249 23474741     Macedo Gde O          0.00 -0.7841640 0.07894525           4  9.359847e-03
# 264 25041182       Tsalikis L          0.00 -0.5555529 0.03318517           3  3.345286e-03
# 266 25222105            Ma WJ          0.00  0.6552411 0.02558360           4  9.889265e-05

# 3) Bootstrapping 
bootstrap_vals = bootstrap_fcn(dataset=hba1c,reliable_values = c(0,1),bootstrap_samples = 1000)
bootstrap_vals

# [[1]]
# value          L        Est          U
# 1 far outlier L -1.2003977 -0.9913397 -0.8436374
# 2     outlier L -0.7096698 -0.5793401 -0.4989443
# 3     outlier U  0.4140373  0.5126489  0.6295490
# 4 far outlier U  0.7616933  0.9237387  1.1170546
# 
# [[2]]
# # A tibble: 4 × 7
# Sigma L_boot_L L_boot L_boot_U U_boot_L U_boot U_boot_U
# <dbl>    <dbl>  <dbl>    <dbl>    <dbl>  <dbl>    <dbl>
#   1     3   -0.670 -0.559   -0.460    0.388  0.515    0.649
# 2     4   -0.896 -0.745   -0.610    0.538  0.701    0.878
# 3     5   -1.13  -0.937   -0.765    0.696  0.893    1.11 
# 4     6   -1.37  -1.14    -0.927    0.856  1.09     1.36 

# 4) Do any of of the identified problematic trials have significant outlier effect sizes 
metaanalysis_effectsizes_lnCVRs <- robu(formula = effect ~  1, data = hba1c[!is.na(hba1c$effect),],
                            studynum = PMID, var.eff.size = var_effect, rho = 0.8,small = TRUE)       
sigma_lvl_diabetes_effectsizes = as.data.frame(cbind(Sigma = c(3,4,5,6),
                                         rbind(robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(3))*2),
                                               robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(4))*2),
                                               robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(5))*2),
                                               robu_predict_interval(metaanalysis_effectsizes_lnCVRs, level = 1-(1-pnorm(6))*2))))
colnames(sigma_lvl_diabetes_effectsizes) = c("Sigma","L","U")
# sigma_lvl_diabetes_effectsizes

# any of detected identified problematic trials have significant outlier effect sizes :
parametric_reliable_df[parametric_reliable_df$sig ==1 & parametric_reliable_df$Outlier != "Not outlier" &
                        abs(parametric_reliable_df$effect/sqrt(parametric_reliable_df$var_effect)) > abs(qnorm(p=0.025)) & 
                         (parametric_reliable_df$effect <  sigma_lvl_diabetes_effectsizes$L[sigma_lvl_diabetes_effectsizes$Sigma==3] |
                            parametric_reliable_df$effect > sigma_lvl_diabetes_effectsizes$U[sigma_lvl_diabetes_effectsizes$Sigma==3]),]

# ---- S4 & S5: baseline characteristics 
library(dplyr); library(officer);library(flextable); library(labelled); library(gtsummary)

set_flextable_defaults(
  font.size = 12, # Adjust font size as needed
  padding = 0,   # Set padding to 0
  line_spacing = 1 # Set line spacing to 1 (or less if desired)
)
list("style_number-arg:big.mark" = "") %>% set_gtsummary_theme()
master <- read_sas("plos.sas7bdat")
master$indexed = ifelse(master$PMID > 0,"PubMed Indexed","Not PubMed Indexed")
master$country_int = ifelse(master$Country %in% c("Egypt","China","India","Iran","Pakistan","Saudi-Arabia","Russia"),
                            "Countries with high retraction rate","Other countries")
# categorical: funding indexed country_int origin_data study_design
# continuous: random_quality blinding_quality ascertainment effect sample_size nGroups duration numb_authors Publication_Year

master  <- master %>% mutate(
    # Study design classification
    study_design = case_when(
      crossoverrho1 == 1 ~ "cross-over",
      crossoverrho1 == 0 ~ "parallel",
      (crossoverrho1 > 0 & crossoverrho1 < 1) | (crossoverrho1 < 0 & crossoverrho1 > -1) ~ "cross-over",
      TRUE ~ NA_character_
    ),
    crossovertrial = case_when(
      study_design == "cross-over" ~ 1,
      study_design == "parallel"  ~ 0,
      TRUE ~ NA_real_
    ),
    # Quality scores (all continuous)
    random_quality   = rowMeans(cbind(r1, r2), na.rm = TRUE),
    blinding_quality = rowMeans(across(b1:b22), na.rm = TRUE),
    ascertainment    = rowMeans(across(a1:a12), na.rm = TRUE),
    baselinePavail = ifelse(is.na(P_onesided),"Unavailable","Available"),
    baseline01temp = ifelse(is.na(P_onesided),NA,ifelse(P_onesided < 0.0005 | P_onesided > 0.9995,1,0)),
    baseline1temp = ifelse(is.na(P_onesided),NA,ifelse(P_onesided < 0.005 | P_onesided > 0.995,1,0)),
    baseline5temp = ifelse(is.na(P_onesided),NA,ifelse(P_onesided < 0.025 | P_onesided > 0.975,1,0))) %>% 
  group_by(PMID,funding,indexed,country_int,origin_data,
            sample_size,nGroups, numb_authors,Publication_Year,
           random_quality,blinding_quality,ascertainment,crossovertrial,
           study_design,baselinePavail) %>% summarize(
             duration=max(floor(month)),effectsize=min(effect,na.rm=T),lncvrsigmamax=max(lncvrsigma),
             baseline01 = max(baseline01temp), baseline1=max(baseline1temp), baseline5 = max(baseline5temp))

master$baseline01 = ifelse(is.na(master$baseline01),NA,ifelse(master$baseline01 == 1,"p<0.001","p>0.001"))
master$baseline1 = ifelse(is.na(master$baseline1),NA,ifelse(master$baseline1 == 1,"p<0.01","p>0.01"))
master$baseline5 = ifelse(is.na(master$baseline5),NA,ifelse(master$baseline5 == 1,"p<0.05","p>0.05"))

# study_design random_quality blinding_quality ascertainment

var_label(master) <- list(funding="Reported funding source",indexed="Indexed in Pubmed",
                          country_int="Country study was conducted in",origin_data="Unpublished data sources",
                          study_design="Cross-over trials", random_quality="Quality score for randomization",
                          blinding_quality="Quality score of blinding",ascertainment="Quality score for ascertainment",
                          effectsize="Max HbA1c effect size",sample_size="Sample size",
                          nGroups="Number of trial arms",duration="Duration of trial",
                          numb_authors="Number of authors",Publication_Year="Publication year",
                          baselinePavail="Data available to calculate\nCarlisle-Fisher-Stouffer p-value",
                          baseline01="Baseline data unbalanced\nor too balanced (p<0.001)",
                          baseline1="Baseline data unbalanced\nor too balanced (p<0.01)",
                          baseline5="Baseline data unbalanced\nor too balanced (p<0.05)")

# 
master_postrandom = as.data.frame(master)
master_postrandom$PMID = as.character(master_postrandom$PMID)
master_postrandom$postrandom = ifelse(is.na(master_postrandom$lncvrsigmamax),"missing","available")
master_postrandom$postrandom = factor(master_postrandom$postrandom,levels=c("missing","available"),
                                      labels=c("missing"="**No standard\ndeviation information**",
                                               "available"="**Standard deviation\ninformation available**"))

master_postrandom |> 
  select(funding,indexed,country_int,origin_data,study_design, 
         random_quality,blinding_quality,ascertainment,effectsize,sample_size,nGroups,
         duration,numb_authors,Publication_Year,
         postrandom,
         baselinePavail,baseline01,baseline1,baseline5) |>  
  tbl_summary(by = postrandom,
              type = list(c(random_quality,blinding_quality,ascertainment,effectsize,sample_size,nGroups,
                            duration,numb_authors,Publication_Year) ~ "continuous2"),
              missing_text = "(Missing)",# missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ({sd})", "{min} - {max}", "{median} ({p25}, {p75})")),
              digits = list(all_continuous() ~ 1,all_categorical() ~ 1)) %>%
  add_p( test =list(all_continuous()  ~ "wilcox.test",all_categorical() ~ "fisher.test" ),
         pvalue_fun = label_style_pvalue(digits = 3))  |>
  bold_labels() |> 
  as_flex_table() %>%autofit() %>%
  set_caption(as_paragraph(as_b("S4 Table - Characteristics of trials by presence or absence of post-randomization standard deviations"))) |> 
  save_as_docx(path="S4_postrandom_sd.docx")

# characteristics of trials with and without problematic baseline data
master_problematicsd = as.data.frame(master)
master_problematicsd$PMID = as.character(master_problematicsd$PMID)

master_problematicsd$problematicsd = ifelse(is.na(master_problematicsd$lncvrsigmamax),"uhoh",
                                            ifelse(master_problematicsd$lncvrsigmamax==1,"Problematic lnCVR",
                                                   "Non-problematic lnCVR"))
master_problematicsd |> filter(problematicsd != "uhoh") |>
  select(funding,indexed,country_int,origin_data,study_design, 
         random_quality,blinding_quality,ascertainment,effectsize,sample_size,nGroups,
         duration,numb_authors,Publication_Year,
         problematicsd,
         baselinePavail,baseline01,baseline1,baseline5) |> 
  tbl_summary(by = problematicsd,
              type = list(c(random_quality,blinding_quality,ascertainment,effectsize,sample_size,nGroups,
                            duration,numb_authors,Publication_Year) ~ "continuous2"),
              missing_text = "(Missing)",# missing = "no",
              statistic = list(all_continuous() ~ c("{mean} ({sd})", "{min} - {max}", "{median} ({p25}, {p75})")),
              digits = list(all_continuous() ~ 1,all_categorical() ~ 1)) %>%
  add_p( test =list(all_continuous()  ~ "wilcox.test",all_categorical() ~ "fisher.test" ),
         pvalue_fun = label_style_pvalue(digits = 2))  |>
  bold_labels() |> 
  as_flex_table() %>% autofit() %>%
  set_caption(as_paragraph(as_b("S5 Table - Characteristics of trials by with and without problematic lnCVR"))) |> 
  save_as_docx(path="S5_problematicsd.docx")








