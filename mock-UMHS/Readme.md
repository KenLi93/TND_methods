---
title: "NC Analysis of Mock U of Michigan Health System Data"
author: "Kendrick Li"
date: '2022-11-11'
output:
  html_document:
    keep_md: yes
---



In this document, we illustrate the proposed NC estimator using a mock University of Michigan Health System data to study SARS-Cov-2 vaccine effectiveness against symptomatic infection.

## Setup and load the dataset

```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 22000)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## system code page: 65001
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] numDeriv_2016.8-1.1 lubridate_1.8.0     MASS_7.3-54        
## [4] dplyr_1.0.8        
## 
## loaded via a namespace (and not attached):
##  [1] knitr_1.37       magrittr_2.0.2   tidyselect_1.1.1 R6_2.5.1        
##  [5] rlang_1.0.1      fastmap_1.1.0    fansi_1.0.2      stringr_1.4.0   
##  [9] tools_4.1.2      xfun_0.34        utf8_1.2.2       cli_3.1.1       
## [13] jquerylib_0.1.4  htmltools_0.5.2  ellipsis_0.3.2   yaml_2.2.2      
## [17] digest_0.6.29    tibble_3.1.6     lifecycle_1.0.1  crayon_1.4.2    
## [21] purrr_0.3.4      sass_0.4.2       vctrs_0.3.8      glue_1.6.1      
## [25] cachem_1.0.6     evaluate_0.14    rmarkdown_2.18   stringi_1.7.6   
## [29] compiler_4.1.2   bslib_0.4.1      pillar_1.7.0     generics_0.1.2  
## [33] jsonlite_1.7.3   pkgconfig_2.0.3
```

```r
tnd_data <- read.csv(file = "mock_tnd_data.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE)
with(tnd_data, table(Vaxname, COVID))
```

```
##          COVID
## Vaxname       0     1
##   Janssen  6249  1220
##   Moderna 15543   887
##   none     9011 10120
##   others    916   208
##   Pfizer   8387   946
```

## Data Preprocessing
1. Categorize age into age groups (younder than 18, 18 to 60, older than 60);
2. Create dummy variables for Charlson Score (0, 1, 2, larger than 3, missing);
3. Create dummy variables for calendar months of the latest COVID test;
4. Create an indicator for Caucasian subjects;
5. Create the negative control outcome (NCO) as a binary indicator of the presence of at least one of the following conditions: arm/leg cellulitis, eye/ear disorder,  gastro-esophageal disease, atopic dermatitis, injuries, and general adult examination visits.

```r
tnd_data_2 <- tnd_data %>%
  mutate(age_le_18 = as.numeric(Age < 18),
         age_geq_18_le_60 = as.numeric(Age >= 18 & Age < 60),
         age_geq_60 = as.numeric(Age >= 60),
         wscore_0 = as.numeric(!is.na(wscore) & wscore == 0),
         wscore_1 = as.numeric(!is.na(wscore) & wscore == 1),
         wscore_2 = as.numeric(!is.na(wscore) & wscore == 2),
         wscore_geq_3 = as.numeric(!is.na(wscore) & wscore >= 3),
         wscore_na = as.numeric(is.na(wscore)),
         april = as.numeric(months == 4),
         may = as.numeric(months == 5),
         june = as.numeric(months == 6),
         july = as.numeric(months == 7),
         aug = as.numeric(months == 8),
         sep = as.numeric(months == 9),
         oct = as.numeric(months == 10),
         nov_dec = as.numeric(months >= 11),
         caucasian = as.numeric(Race == "Caucasian"),
         nco = as.numeric(EE + CLT + gfd + dermatitis + injuries + GE_NCO > 0))
```

## Detecting the confounding bias
1. Adjusted logistic regression of the NCO on vaccination and NCE 

```r
glm_NCO <- glm(nco ~ vaccinated + IMM + age_geq_18_le_60 + age_geq_60 + Gender + caucasian + wscore_geq_3 + april + may + june + july + aug + sep + oct,
               family = binomial,
               data = tnd_data_2)
summary(glm_NCO)
```

```
## 
## Call:
## glm(formula = nco ~ vaccinated + IMM + age_geq_18_le_60 + age_geq_60 + 
##     Gender + caucasian + wscore_geq_3 + april + may + june + 
##     july + aug + sep + oct, family = binomial, data = tnd_data_2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.2804  -1.1635  -0.9918   1.1786   1.4004  
## 
## Coefficients:
##                   Estimate Std. Error z value Pr(>|z|)    
## (Intercept)      -0.389076   0.042133  -9.234  < 2e-16 ***
## vaccinated        0.398995   0.018690  21.348  < 2e-16 ***
## IMM               0.137311   0.017961   7.645 2.09e-14 ***
## age_geq_18_le_60 -0.013433   0.034587  -0.388   0.6977    
## age_geq_60       -0.043408   0.035019  -1.240   0.2151    
## Gender            0.013421   0.017487   0.767   0.4428    
## caucasian        -0.012742   0.020458  -0.623   0.5334    
## wscore_geq_3      0.021788   0.028047   0.777   0.4373    
## april             0.026115   0.032455   0.805   0.4210    
## may               0.069163   0.032256   2.144   0.0320 *  
## june              0.002863   0.032575   0.088   0.9300    
## july              0.043713   0.032145   1.360   0.1739    
## aug              -0.012345   0.031970  -0.386   0.6994    
## sep              -0.065180   0.030956  -2.106   0.0352 *  
## oct              -0.037403   0.030872  -1.212   0.2257    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 74016  on 53486  degrees of freedom
## Residual deviance: 73405  on 53472  degrees of freedom
## AIC: 73435
## 
## Number of Fisher Scoring iterations: 4
```
The p-values for the coefficients of SARS-Cov-2 vaccination and NCE (prior immunization) are <2e-16 and 2.09e-14, respectively, indicating the presence of confounding bias.

2. Adjusted logistic regression of COVID-19 on vaccination and NCE 

```r
glm_NCE <- glm(COVID ~ vaccinated + IMM + age_geq_18_le_60 + age_geq_60 + Gender + caucasian + wscore_geq_3 + april + may + june + july + aug + sep + oct,
               family = binomial,
               data = tnd_data_2)
summary(glm_NCE)
```

```
## 
## Call:
## glm(formula = COVID ~ vaccinated + IMM + age_geq_18_le_60 + age_geq_60 + 
##     Gender + caucasian + wscore_geq_3 + april + may + june + 
##     july + aug + sep + oct, family = binomial, data = tnd_data_2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.6001  -0.5397  -0.4135   0.8072   2.5659  
## 
## Coefficients:
##                   Estimate Std. Error z value Pr(>|z|)    
## (Intercept)       0.365377   0.054087   6.755 1.43e-11 ***
## vaccinated       -2.399133   0.024593 -97.553  < 2e-16 ***
## IMM              -0.349143   0.024743 -14.111  < 2e-16 ***
## age_geq_18_le_60  0.027552   0.045558   0.605    0.545    
## age_geq_60        0.116805   0.046504   2.512    0.012 *  
## Gender            0.408280   0.023499  17.375  < 2e-16 ***
## caucasian         0.008773   0.027048   0.324    0.746    
## wscore_geq_3      0.055110   0.037503   1.469    0.142    
## april            -0.871125   0.044728 -19.476  < 2e-16 ***
## may              -0.723347   0.043605 -16.589  < 2e-16 ***
## june             -0.822736   0.044666 -18.420  < 2e-16 ***
## july             -0.827355   0.044262 -18.692  < 2e-16 ***
## aug              -0.813729   0.043895 -18.538  < 2e-16 ***
## sep              -0.045710   0.039260  -1.164    0.244    
## oct              -0.061323   0.039291  -1.561    0.119    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 60176  on 53486  degrees of freedom
## Residual deviance: 46318  on 53472  degrees of freedom
## AIC: 46348
## 
## Number of Fisher Scoring iterations: 5
```
The p-value for the coefficient NCE (prior immunization) is <2e-16, also indicating the presence of confounding bias.

## NC estimator of vaccine effectiveness
1. Pfizer-BioNTech vaccine

```r
pfizer_data <- tnd_data_2 %>% filter(Vaxname %in% c("Pfizer", "none"))

pfizer_iptw <- tnd_iptw(data = list(A = pfizer_data$vaccinated,
                                    Y = pfizer_data$COVID,
                                    X = pfizer_data[, c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3", "april", "may", "june", "july", "aug", "sep", "oct")]))

pfizer_nc <- nc_covid_ve(data = pfizer_data,
                         nce = "IMM",
                         nco = "nco",
                         x = c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3", "april", "may", "june", "july", "aug", "sep", "oct"))
```
The estimated VE and their 95\% confidence intervals are:

* Logistic regression: 0.907 (0.9, 0.914)
* IPTW: 0.714 (0.695, 0.731)
* NC estimator: 0.86 (0.832, 0.883)


2. Moderna vaccine

```r
moderna_data <- tnd_data_2 %>% filter(Vaxname %in% c("Moderna", "none"))

moderna_iptw <- tnd_iptw(data = list(A = moderna_data$vaccinated,
                                    Y = moderna_data$COVID,
                                    X = moderna_data[, c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3", "april", "may", "june", "july", "aug", "sep", "oct")]))

moderna_nc <- nc_covid_ve(data = moderna_data,
                         nce = "IMM",
                         nco = "nco",
                         x = c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3", "april", "may", "june", "july", "aug", "sep", "oct"))
```
The estimated VE and their 95\% confidence intervals are:

* Logistic regression: 0.953 (0.949, 0.956)
* IPTW: 0.807 (0.794, 0.819)
* NC estimator: 0.892 (0.833, 0.931)


3. Johnson \& Jognson's Janssen vaccine

```r
janssen_data <- tnd_data_2 %>% filter(Vaxname %in% c("Janssen", "none"))

janssen_iptw <- tnd_iptw(data = list(A = janssen_data$vaccinated,
                                    Y = janssen_data$COVID,
                                    X = janssen_data[, c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3", "april", "may", "june", "july", "aug", "sep", "oct")]))

janssen_nc <- nc_covid_ve(data = janssen_data,
                         nce = "IMM",
                         nco = "nco",
                         x = c("age_geq_18_le_60", "age_geq_60", "Gender", "caucasian", "wscore_geq_3", "april", "may", "june", "july", "aug", "sep", "oct"))
```
The estimated VE and their 95\% confidence intervals are:

* Logistic regression: 0.839 (0.828, 0.85)
* IPTW: 0.57 (0.546, 0.592)
* NC estimator: 0.74 (0.669, 0.795)
