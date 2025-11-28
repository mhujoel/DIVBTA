# DIVBTA for Boldt et al. (2009)

To assess whether the proposed methodology would identify the study of Boldt et al. as potentially problematic we would need to both:
* identify the DIVBTA from Boldt et al. as being statistically *significant* 
* identify the DIVBTA from Boldt et al. as having an *outlying* value 

The mean, standard deviation, and sample sizes for the two arms can be obtained from the abstract of Boldt et al. (DOI: 10.1213/ANE.0b013e3181b5a24b): "Base excess after surgery was lower in the albumin-based priming group than in the balanced HES priming group (-5.9 +/- 1.2 mmol/L vs +0.2 +/- 0.2 mmol/L, P = 0.0003)." 

All analyses are provided in `Boldt_DIVBTA.R`; but some highlights below: 

**Significance of DIVBTA from Boldt et al** 

If we test the significance of the lnVR from Boldt et al:

```
xT=-5.9;sT=1.2;nT=25;xC=0.2;sC=0.2;nC=25
# F-test
ratio_var = sT^2/sC^2
(2 * ifelse(ratio_var < 1,
  pf(ratio_var, df1=nT-1, df2=nC-1, lower.tail = TRUE),
  pf(ratio_var, df1=nT-1, df2=nC-1, lower.tail = FALSE)))

[1] 3.111791e-13

lnVR = log(sT/sC)+1/(2*(nT-1))-1/(2*(nC-1))
sd_lnVR = sqrt(1/(2*(nT-1))+1/(2*(nC-1)))
# T-test:
(2*exp(pt(abs(lnVR/sd_lnVR),df=nT+nC-2,lower.tail = F,log.p = T)))

[1] 1.509264e-11

# Assumes normality:
(2*exp(pnorm(abs(lnVR/sd_lnVR),lower.tail = F,log.p = T)))

[1] 1.667139e-18
```

The DIVBTA from Boldt et al. **is** statistically significant (p=3.1e-13). 

**Is the DIVBTA from Boldt et al an outlier?** 

We used the data from Figure 2 of Curran JD et al. Comparison of Balanced Crystalloid Solutions: A Systematic Review and Meta-Analysis of Randomized Controlled Trials. *Crit Care Explor.* 2021;3(5):e0398 which had base excess values from randomized controlled trials (provided in `base.csv`). We used the presented 10 values with one modification:
* Benoit 2016: the standard deviation for the comparator was changed from 0.7 to 1.26; while the SD for the Plasmalyte for Benoit was correctly estimated from the quartile range as 0.96 (SD ~ IQR/1.35 = 1.3/1.35 = 0.96) the SD for the comparator for Benoit, when derived the same way, should be 1.26 not 0.07 (1.7/1.35 = 1.26). 

The lnVR for Boldt et al. can be computed to be:
```
> xT=-5.9;sT=1.2;nT=25;xC=0.2;sC=0.2;nC=25
> (lnVR = log(sT/sC)+1/(2*(nT-1))-1/(2*(nC-1)))
[1] 1.791759
```

Tukey Fences for the lnVR across the 10 studies can be found to be: 
```
   extreme_L     mild_L    mild_U extreme_U
1 -0.8208672 -0.5381081 0.2159162 0.4986753
```
accordingly, the lnVR for Boldt et al. is an *extreme* outlier using a non-parametric approach.

Using a parametric approach and computing prediction intervals for the lnVR across the 10 studies:
```
     Sigma level          L         U
[1,]           3 -0.8837272 0.6239395
[2,]           4 -1.4686589 1.2088712
[3,]           5 -2.6416245 2.3818368
[4,]           6 -5.3014579 5.0416701
```
accordingly, the lnVR for Boldt et al. is an *extreme* outlier using a parametric approach (4-sigma event).



