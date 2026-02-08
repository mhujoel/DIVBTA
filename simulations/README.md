# Simulation Studies

R libraries required: sn; dplyr; fitdistrplus; ggplot2; doParallel; foreach; grid; magick; pdftools; haven; parallel; patchwork.

# How we define 3-sigma outlier:

We use the data from Section 3 to estimate bounds for a 3-sigma outlier lnCVR for our simulation study:

```
library(robumeta); 
diabetes_trials=read.csv("github_LnCVR.csv")
# restrict to trials we are more confident in:
diabetes_trials_confident = diabetes_trials[diabetes_trials$P_onesided > 0.025 &
                                              diabetes_trials$P_onesided < 0.975 &
                                              !is.na(diabetes_trials$P_onesided), ]
# perform our meta-analysis using robumeta:
fit <- robu(formula = LnCVR ~ 1,data=diabetes_trials_confident , 
                     studynum = PMID, var.eff.size=var_LnCVR)

level = 1-(1-pnorm(3))*2 # 3-sigma event
# in simulations we assume no effect thus M*=0 

# For df we can use N-2 as suggested by Introduction to Meta-Analysis (2nd edition)
# Chapter 17 page 122 or fit$dfs provided by robumeta

# we can test both df and the 3-sigma bounds they result in:
# df we can use N-2 as suggested by Introduction to Meta-Analysis (2nd edition)
# Chapter 17 page 122
(df_val = fit$N - 2) ; p_crit <- (1 - level) / 2
# [1] 173 <- 173 degrees of freedom 
0 + c(-1, 1) *qt(p = p_crit, df = df_val, lower.tail = FALSE) *
           sqrt(as.numeric(fit$mod_info$tau.sq) + as.numeric(fit$reg_table[3])^2)
# [1] -0.5001689  0.5001689

# fit$dfs provided by robumeta
(df_val = fit$dfs) ; p_crit <- (1 - level) / 2
# [1] 143.3904 <- ~143 degrees of freedom
0 + c(-1, 1) *qt(p = p_crit, df = df_val, lower.tail = FALSE) *
  sqrt(as.numeric(fit$mod_info$tau.sq) + as.numeric(fit$reg_table[3])^2)
# [1] -0.5016838  0.5016838

# based on the above we use -0.5 and 0.5 to define 3-sigma outliers in our simulation studies
```



# Heterogeneous treatment effects
Heterogeneous treatment effects (HTE) are defined as a variation in treatment effect in the intervention arm of a trial which is not random but systematically related to observable or unobservable patient characteristics such as age, genetic markers, or baseline biomarker levels. The impact of HTE on DiVBTAs was assessed by modeling two parameters: (1) the proportion of trial participants in the treatment group who respond to the intervention (this proportion was varied between 0% and 50%), and (2) the magnitude of the treatment effect in the proportion of trial participants responding to treatment (which was varied from a 0.2% improvement to a 1.4% improvement in HbA1c level). 

# Missing-non-at-random dropout
Missing-non-at-random dropout is a type of missing data mechanism in which the probability that a participant drops out (i.e., their outcome data are missing) depends on the unobserved (missing) outcome values themselves, even after accounting for all observed data. For instance, participants on a downhill clinical trajectory may drop out of trials when realizing they are in the control group. To model dropout, the worst responders were deleted from the control group.  This type of dropout can be considered extreme, unlikely to occur in real-world clinical trial settings. The proportion of worst responders dropping from the control group was varied between 0% to 50%.  

# Duplication of data for good responders
The aim of including these simulations was to assess whether statistically significant DiVBTAs outliers could offer a sensitive statistic to detect data manipulation in the form of data duplication.  The impact of duplicating up to 20% of the best responders in the intervention group on LnCVR was investigated in small and large samples. 
