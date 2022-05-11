# Single-Cell Differential Expression Analysis (sc-DEA) with Telomere Length (TL)
We provide three main scripts for the scDEA with TL: 

**1. [dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/dea_MAST_glmer_TL.R):** to perform the sc-DEA with TL per cell type at the high (approach I) or low (approach II) resolution level. 

**2. [dea_MAST_statistics.R](Projects/Telomere_analysis/Single_cell_analysis/dea_MAST_statistics.R):** to summarize the sc-DEA with TL results at the cell type resolution level (approach I or approach II). 

**3. [Downstream_SC_DE.R](Projects/Telomere_analysis/Single_cell_analysis/Downstream_SC_DE.R):** to peform the downstream analysis on the sc-DEA with TL from the approach II.

*Of note*: The functions used in *dea_MAST_glmer_TL.R* and *dea_MAST_statistics.R* are defined in an [accessory script](Projects/Telomere_analysis/Single_cell_analysis/scripts/accessory_functions.R).

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es) or Sergio Andreu-Sánchez (sergioandreusanchez@gmail.com). 


## Required Software
* **R** >=4.0.0 version: You need to install the packages loaded in the [accessory script](scripts/accessory_functions.R).

*Of note*: To install the [MAST R package](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5), please use the following command to ensure you install the MAST version >= 1.16.0: 

``` 
devtools::install_github("RGLab/MAST")
packageVersion("MAST") #>=1.16.0
```

## Required Input
This section explains the input data and it’s structure to run the three main scripts. To follow the explanations in the **Required Input** section, you can clone this repository and change your current working directory.   

```
git clone \
  --depth 3  \
  --filter=blob:none  \
  --sparse \
  https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome \
;
cd Groningen-Microbiome
git sparse-checkout set Projects/Telomere_analysis/Single_cell_analysis
cd Projects/Telomere_analysis/Single_cell_analysis/
```

*Of note*: `git clone --filter` from git 2.19 is used to clone only objects in the [Projects/Telomere_analysis/Single_cell_analysis](Projects/Telomere_analysis/Single_cell_analysis) directory. 

### Test Data
We have provided a **cell-type-specific example dataset** for each of the **two approaches**: CD8T memory cells for approach I (n =  11,071) and CD8T cells for approach II (n =  13,246), together with some other inputs required to run the sc-DEA with TL [dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/dea_MAST_glmer_TL.R). These files are hosted in the [**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/) in the [MOLGENIS cloud](https://www.molgenis.org/).

Here is the structure of the **input directory (path/to/in_dir)** for [dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/dea_MAST_glmer_TL.R). You can find this same structure in the input directory example in [**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/):

**input directory (path/to/in_dir)**    
├── approach_I  
│   └── CD8T_memory.rds  
├── approach_II  
│   └── CD8Tcells.rds  
├── covariates.approach_II.tab  
├── covariates.approach_I.tab  
├── top10genes.approach_II.tab  
├── top10genes.approach_I.tab  
├── top3genes.approach_II.tab  
└── top3genes.approach_I.tab  

*Of note*: Some of the files are strictly required, but some of them are optional:
**Required**: 
* You should have a directory for each of the cell type classification levels (e.g., approach_I and approach_II directories) containing the seurat object for a specific cell type (e.g., approach_I > CD8T_memory.rds, and approach_II > CD8Tcells.rds).
* You should have a tab separated file with the covariates considered in the scDEA model (e.g., covariates.approach_I.tab and covariates.approach_II.tab). The covariates names has to be the same in the metadata slot of the seurat object. The details of this file will be explained in the next section.

**Optional**: 
* Since the scDEA with TL is peformed at the gene-wise level across all the expressed genes, it can take a lot of time to run depending on the number of cells in you cell-type-specific seurat object. Thus, we provide several tab files to try [dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/dea_MAST_glmer_TL.R) only with a subset of genes (e.g., top3genes.approach_I.tab, top3genes.approach_II.tab, top10genes.approach_I.tab, and top10genes.approach_II.tab). The details of this file will be explained in the next section.

### Required Data
**input directory (path/to/in_dir)**    
├── approach_I  
│   └── CD8T_memory.rds  
├── approach_II  
│   └── CD8Tcells.rds  
├── covariates.approach_II.tab  
├── covariates.approach_I.tab    

#### QC-MAD combinations ([qc_mad.tab](qc_mad.tab))
A tsv file that has in the:
* 1st column: QC metrics. By default, number of UMIs (*nCount_RNA*) and % of mitochondrial genes (*percent.mt*).
* 2nd column: Upper or lower threshold. By default, lower for *nCount_RNA* and upper for *percent.mt*.
* 3rd and 4rd columns: minimum and maximum MADs. By default, *minimum*=1 and *maximum*=5.

*Of note*:
* Tab separated
* It is assumed that the QC metrics are calculated in the seurat object as a result from WG1 pipeline, and thus, they are columns of the metadata slot of the seurat object.
* This file must have this header. 
* The QC-MAD combinations file provided for the test dataset is the [qc_mad.tab](/qc_mad.tab) file:

| QC_metric  | bound | MAD_min  | MAD_max |
| ------------- | ------------- | ------------- | ------------- |
| nCount_RNA  | lower  | 1  | 5 |
| percent.mt  | upper  | 1  | 5 |

#### Azimuth l1-l2 pairing file ([azimuth_l1_l2.csv](/azimuth_l1_l2.csv))
A csv file that has in the:
* 1st column: Azimuth's level 1 cell type classification (L1).
* 2nd column: Azimuth's level 2 cell type classification (L2).

### Optional Data
**input directory (path/to/in_dir)**    
├── top10genes.approach_II.tab  
├── top10genes.approach_I.tab  
├── top3genes.approach_II.tab  
└── top3genes.approach_I.tab

#### Metadata variables ([metadata_variables.tab](/metadata_variables.tab))
A tsv file that has in the:
* 1st column: Metadata variable name. 
* 2nd column: Metadata variable type. 
* 3rd and 4rd columns: minimum and maximum MADs. By default, *minimum*=1 and *maximum*=5.

*Of note*:
* Tab separated.
* It is assumed that the metadata variable names are columns of the metadata file or metadata slot of the seurat object.
* This file must have this header.
* By default, the QC statistics will be summarized at the whole dataset. You can choose to summarize them by metadata variable.
* The metadata variables file provided for the test dataset is the [metadata_variables.tab](/metadata_variables.tab) file: 

| md_var  | type |  
| ------------- | ------------- |  
| Pool  | donor  |  
| Assignment  | donor  |  
| predicted.celltype.l2  | cell  |  
| scpred_prediction  | cell  |  
| predicted.celltype.l1  | cell  |  

## Running the [add-on script](QC_statistics.R)
*Of note*: The functions called in the [add-on script](QC_statistics.R) are defined in an [external script](scripts/QC_functions.R).

**1.** If you have not done it yet, the first step would be to clone this repository and change your current working directory.    
```
git clone https://github.com/aidarripoll/wg1-qc_filtering.git  
cd wg1-qc_filtering
```

**2.** Set common environmental variables:  
```
dataset=wg2_onek1k_subset  
input_directory=wg2-cell_type_classification
output_directory=QC_statistics_examples
```

**3.** Running the add-on script with different parameters:  

3.1. Summarize the QC statistics at the dataset level after:

* Calculating the QC statistics at the dataset level:
```
Rscript QC_statistics.R --dataset $dataset --in_dir $input_directory --out_dir $output_directory
```

* Calculating the QC statistics at the batch metadata variable (i.e., Pool) level:
```
batch_variable=Pool  
Rscript QC_statistics.R --dataset $dataset --level $batch_variable --in_dir $input_directory --out_dir $output_directory 
```

The output for each of the parameters settings is the summary statistics (*tag.rds*):
QC_statistics_examples  
└── wg2_onek1k_subset  
    └── nCount_RNA_lower_1_5.percent.mt_upper_1_5  
        └── **by_dataset**   
            ├── dataset  
            │   └── **tag.rds**  
            └── Pool  
                └── **tag.rds**  
                
# Example outputs

To reproduce the supplementary tables for the example cell types: CD8T memory cells **Table_S7.1 (approach I)** and CD8T cells **Table_S7.2 (approach II)**, you should run the following commands:

* CD8T memory cells **Table_S7.1 (approach I):**

* CD8T cells **Table_S7.2 (approach II)**:



We have provided the two output directories for the test data *(wg2_onek1k_subset)* in a tar.gz format:

* [QC_statistics_examples.tar.gz](QC_statistics_examples.tar.gz): Outputs from running the commands in **3.1, 3.2 and 3.3** as part of **Running the add-on script** section.

* [QC_statistics_examples.files.tar.gz](QC_statistics_examples.files.tar.gz): Outputs from running the commands in **3.1 and 3.2** as part of **Running the add-on script** section, and also running the **alternative** option in the **Discussion in the QC Threshold Selection Committee** section.


You can decompress them by:
```
tar -xvf QC_statistics_examples.tar.gz
tar -xvf QC_statistics_examples.files.tar.gz
```

