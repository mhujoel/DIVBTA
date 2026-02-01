# Simulation Studies

R libraries required: sn; dplyr; fitdistrplus; ggplot2; doParallel; foreach; grid; magick; pdftools; haven; parallel; patchwork.

# Heterogeneous treatment effects
Heterogeneous treatment effects (HTE) are defined as a variation in treatment effect in the intervention arm of a trial which is not random but systematically related to observable or unobservable patient characteristics such as age, genetic markers, or baseline biomarker levels. The impact of HTE on DiVBTAs was assessed by modeling two parameters: (1) the proportion of trial participants in the treatment group who respond to the intervention (this proportion was varied between 0% and 50%), and (2) the magnitude of the treatment effect in the proportion of trial participants responding to treatment (which was varied from a 0.2% improvement to a 1.4% improvement in HbA1c level). 

# Missing-non-at-random dropout
Missing-non-at-random dropout is a type of missing data mechanism in which the probability that a participant drops out (i.e., their outcome data are missing) depends on the unobserved (missing) outcome values themselves, even after accounting for all observed data. For instance, participants on a downhill clinical trajectory may drop out of trials when realizing they are in the control group. To model dropout, the worst responders were deleted from the control group.  This type of dropout can be considered extreme, unlikely to occur in real-world clinical trial settings. The proportion of worst responders dropping from the control group was varied between 0% to 50%.  

# Duplication of data for good responders
The aim of including these simulations was to assess whether statistically significant DiVBTAs outliers could offer a sensitive statistic to detect data manipulation in the form of data duplication.  The impact of duplicating up to 20% of the best responders in the intervention group on LnCVR was investigated in small and large samples. 
