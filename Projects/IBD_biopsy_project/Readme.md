# Mucosal host–microbe interactions associate with clinical phenotypes in inflammatory bowel disease

This project is aiming to identify intestinal gene expression and microbiota change in repsect to the heterogeneity of IBD. We also want to highlight that considering the host-microbiota interactions would reveal more novel insights into the disease putative mechanisms compared solely with transcriptomic or microbial data.

### Raw data processing
Pipeline folder

### Host intestinal transcriptomic and microbial characterization

1) ProjectCode.RNAseq.R
2) ProjectCode.16S.R
3) HALLA.R (https://github.com/biobakery/halla_legacy)
4) PairedSample.Analysis.R

### Host-microbiota interactions

1) SparseCCA.R
2) Network.analysis.R
3) Dysbiosis.R
4) Deconvolution.R

***
All the analysis have been adjusted for potential confounders if neccessary, including age, sex, BMI, tissue location and inflammation, medication use(aminosalicylates, thiopurines and steroids),sample batch and surgery (for 16S data only). The issue multiple samples from one patient is taken into account by introducing a random effect (1|ID) into the generalized linear mixed model.
***

### Data deposit
All the data is at EGA (https://ega-archive.org/studies/EGAS00001002702).
### Analysis support
Please contact dhu.sxhu@hotmail.com;a.r.bourgonje@umcg.nl
