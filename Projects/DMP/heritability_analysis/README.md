## DMP heritability analysis, POLY models and workflow

The folder contains scripts used to execute heritability analysis and POLY model files

### Links to software used in heritability analysis:

- POLY: https://people.virginia.edu/~wc9c/poly/input.html
- pedtools: https://csg.sph.umich.edu/abecasis/pedstats/download/index.html

### Files:

- POLY_do_heritability.sh : bash script for performing heritability analysis, includes cleaning and trimming of pedigree files, calculation of pedigree statistics and summary statistics for traits, and heritability analysis using additive model
- __model_additive.model : POLY heritabiltiy model
- ___polyrun_dag3_clrnorm.txt : POLY output
- prepHeritabilityPlots_v3 : R script for heritability plots, requires trait description file (./data), heritability results in ./data and scripts in ../r_scripts_library
- prepBcDissimilarityPlots_v3.R : R script for plots of results of pair-wise Bray-Curtis dissimilarity analysis
- ./data : input data for prepHeritabilityPlots_v3.R
- ./Plots : output from prepHeritabilityPlots_v3.R and prepBcDissimilarityPlots_v3.R
