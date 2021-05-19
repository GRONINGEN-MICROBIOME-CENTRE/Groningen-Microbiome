## DMP heritability analysis, heritability analysis (V2)

The folder contains R scripts used to execute heritability analysis for DMP project, mock data for testing of codes, and expected results from running codes on the mock data.

### Notes:

- mock_data folder contains *mock* data, this data is intended to be used only for testing of codes. It was designed to produce results similar to real data, but it is *not* real DMP data. 
- real microbiome data used for DMP project can be obtained from EGA study EGAS00001005027, real kinship matrices, family information and cohousing data can be obtained from Lifelines Biobank
- example results generated from mock data are provided (zipped) in mock_data_heritability_results. These results are expected output from DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R
- results permutation runs are NOT provided due to github space constrains, these can be generated using DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R (see Instructions)
- example consolidated results produced from mock data are provided in root folder. These are expected results from running DMP_heritability_v10_mockdata_collectdata.R on mock data heritability models
- example plots generated from mock data results are provided in Plots folder. These are expected results from running DMP_heritability_v10_mockdata_plotresults.R
- codes were developed & tested on R 3.6.0 (see R_session_info.txt for list of used packages)

### Instructions for code testing:

- download the codes and data
- unzip the mock data (mock_data/*.zip)
- Set working folders to correct locations in the scripts (DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R [line 70]; DMP_heritability_v10_mockdata_plotresults.R [line 42]; DMP_heritability_v10_mockdata_collectdata.R [line 15]
- run DMP_heritability_v10_mockdata_taxa.R to generate heritability model for taxa and DMP_heritability_v10_mockdata_pwy.R for pathways
- run DMP_heritability_v10_mockdata_collectdata.R to consolidate data and calculate study-wide FDRs from permutation runs
- run DMP_heritability_v10_mockdata_plotresults.R to generate heritability plots
- NOTE: DMP_heritability_v10_mockdata_plotresults.R can be run directly on the provided consolidated data, example input data is provided
- NOTE2: DMP_heritability_v10_mockdata_collectdata.R requires premutation runs produced by MP_heritability_v10_mockdata_taxa.R and DMP_heritability_v10_mockdata_pwy.R 

### Files:

- README.me : this readme file
- bacpaths.txt : list of bacterial pathways used in DMP project (required by DMP_heritability_v10_mockdata_collectdata.R script)
- DMP_heritability_v10_mockdata_taxa.R : script for heritability analysis of microbiome taxa, requires unzipped mock data
- DMP_heritability_v10_mockdata_pwys.R : script for heritability analysis of microbiome pathways, requires unzipped mock data
- DMP_heritability_v10_mockdata_collectdata.R : script for collecting results and study-wide FDR analysis using permutation runs. Requires results of DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwys.R
- DMP_heritability_v10_mockdata_plotresults.R : script for plotting heritability results. Requires results of DMP_heritability_v10_mockdata_collectdata.R
- results_mockdata_taxa.csv : results of heritability analysis of taxa
- results_mockdata_pwys.csv : results of heritability analysis of pathways
- results_mockdata_withFDRs_and_CIs_taxa.csv : results of heritability analysis of taxa with FDRs calculated from permutation runs
- results_mockdata_withFDRs_and_CIs_pwys.csv : results of heritability analysis of pathways with FDRs calculated from permutation runs
- R_session_info.txt: list of R packages and versions used in development & testing
- ./mock_data : mock data for code testing (zipped)
- ./mock_data_heritability_results : results of heritability analysis of mock data (zipped)
- ./mock_data_permutation_runs : example results of permutation runs of one taxon and one pathway
- ./Plots : heritability plots generated from heritability analysis of mock data