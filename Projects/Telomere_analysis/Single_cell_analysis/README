# Single-Cell Differential Expression Analysis (sc-DEA) with Telomere Length
We provide three main scripts for the scDEA with telomere length: 

**1**. dea_MAST_glmer_TL.R: to perform the sc-DEA with telomere length per cell type at the high (approach I) or low (approach II) resolution level. 

**2**. dea_MAST_statistics.R: to summarize the sc-DEA with telomere length results at the cell type resolution level (approach I or approach II). 

**3**. Downstream_SC_DE.R: to peform the downstream analysis on the sc-DEA with telomere length from the approach II.

*Of note*: The functions used in *dea_MAST_glmer_TL.R* and *dea_MAST_statistics.R* are defined in an [accessory script](scripts/accessory_functions.R).

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es) or Sergio Andreu-Sánchez (sergioandreusanchez@gmail.com). 


## Required Software
* **R** >=4.0.0 version: You need to install the packages loaded in the [accessory script](scripts/accessory_functions.R).

*Of note*: To install the [MAST R package](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5), please use the following command to ensure the MAST version is >= 1.16.0: 

``` 
devtools::install_github("RGLab/MAST")
packageVersion("MAST") #>=1.16.0
```

## Required Input
This section explains the input data and it’s structure to run the three main scripts.

*Of note*: To follow better the explanations in the **Required Input** section, you can clone this repository and change your current working directory.   

```
git clone https://github.com/aidarripoll/wg1-qc_filtering.git  
cd wg1-qc_filtering
```

### Test Data
We have provided a **test dataset** *(wg2_onek1k_subset)* that contains one pool of a 10x run from the [**OneK1K** dataset](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02293-3). Notice that it is a significantly down-sized and sub-sampled version of the whole dataset. In this test dataset, the total number of cells is 1,207 from 13 donors.

Here is the structure of the [input directory for the test dataset](/wg2-cell_type_classification/wg2_onek1k_subset/). This input directory (*/wg2-cell_type_classification/wg2_onek1k_subset/*) should have the same structure as the WG2 pipeline output directory. We will need only the files in the [step4_reduce](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/) directory:

**wg2-cell_type_classification**    
└── wg2_onek1k_subset  
    ├── cell_classification.sif  
    ├── map_hierscpred.R  
    ├── schier_workaroung.sh  
    ├── step1_split  
    │   └── OneK1K-test_dataset.RDS  
    ├── step2_azimuth  
    │   ├── OneK1K-test_dataset.RDS  
    │   ├── OneK1K-test_dataset_ref_spca.png  
    │   └── OneK1K-test_dataset_ref_umap.png  
    ├── step3_hierscpred  
    │   └── OneK1K-test_dataset.RDS  
    ├── **step4_reduce**   
    │   ├── **metadata.reduced_data.RDS**    
    │   └── **reduced_data.RDS**    
    └── step5_compare  
        ├── comparison_contingency_table.tsv  
        ├── comparison_heatmap_counts.pdf  
        └── comparison_heatmap_prop.pdf  
        
The main input for the [add-on script](QC_statistics.R) is the metadata slot ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)) of the seurat object provided by WG2 pipeline ([reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/reduced_data.RDS)). The WG2 pipeline is peforming the cell type classification of the non-QC filtered singlets predicted by WG1 pipeline.

* **Recommended:** We recommend you to use the WG2 seurat object's metadata slot ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)). It will improve the running time and memory of the script. 

* **Alternative:** You can also use the WG2 seurat object ([reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/reduced_data.RDS)). However, it will slow down the running time and memory of the script as we will need to read the full seurat object which can be very large depending on the number of cells (e.g., ~77K cells, 8.9G). 

*Of note*: 
* At this moment, the WG2 pipeline is not providing the ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)) yet. Although you can run the [add-on script](QC_statistics.R) using the whole seurat object, we encourage you to save the metadata slot with the name *metadata.reduced_data.RDS* in the step4_reduce/ directory provided by WG2 pipeline before running the [add-on script](QC_statistics.R) to improve the running time and memory of the script. 

* In case your dataset contains V2 and V3 chemistries, you should create different metadata or seurat objects files in order to run this [add-on script](QC_statistics.R) separately. If this had been the case of this test dataset, you would have ended up with two different datasets (e.g., wg2_onek1k_subset.V2 and wg2_onek1k_subset.V3), meaning that the [add-on script](QC_statistics.R) would have been run separately in each of these two datasets.
