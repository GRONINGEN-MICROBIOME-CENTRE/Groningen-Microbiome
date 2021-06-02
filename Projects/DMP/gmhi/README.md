## GMHI analysis

This folder contains scripts used to calculate Gut Microbiome Health Index (GMHI) on the DMP data. GMHI was introduced in study by Gupta et al. (Nature Comm. 2020; https://www.nature.com/articles/s41467-020-18476-8 & https://github.com/jaeyunsung/GMHI_2020)

### Files:

- gmhi_schematic.png : schematic of the GMHI calculation
- gmhi.R : implementation of GMHI algorithm
- run_gmhi2.R : script that runs gmhi algorithm on DMP data

# 1 We ran the GMHI classifier to see how well we can predict disease status ("healthy" vs "unhealthy") based on our samples. We used a 90-10 training-testing split to generate a discovery (training data) and validation (testing data).  
# 2 We ran the GMHI classifier on all our data using the optimized parameters found by Gupta et al. (2020).        

Breifly, the GMHI classifier is first run the discovery cohort (training data) to find the parameters which optimizes the balanced accuracy, i.e. the prevalence fold change (theta_f) and the prevalence difference (theta_d). Then the validation cohort (testing data) is then run through the classifier, and the balanced accuracy corresponding to the optimized parameters (theta_f and theta_d from the training step) is extracted.   
 
