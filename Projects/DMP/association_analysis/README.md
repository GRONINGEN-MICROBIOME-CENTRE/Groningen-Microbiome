## DMP microbiome variance analysis

The folder contains scripts used to calculate microbiome variance explains by phenotypes
and groups of phentypes

### Files:

- DMP_do_univar_adonis_taxa.R : script for quantification of microbiome composition variance explained by one phenotype
- DMP_do_univar_adonis_PWYs.R : script for quantification of microbiome functionality variance explained by one phenotype
- DMP_do_multivar_adonis.R : script for quantification of variance of microbiome functionality or composition explained by a group of phenotypes
- DMP_collect_adonis_results.R : script for collection of results produced by DMP_do_univar_adonis_taxa or DMP_do_univar_adonis_PWYs results
- DMP_check_adonis_results.R : script for identification of failed adonis jobs
- DMP_variance_explained_plots.R : script for plotting results of multivariate adonis analysis
- ./group_adonis/ : folder with input files for multivariate adonis plots and output of the script
