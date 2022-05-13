# Single-Cell Differential Expression Analysis (sc-DEA) with Telomere Length (TL)
We provide three main scripts for the scDEA with TL: 

**1. [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R):** to perform the sc-DEA with TL per cell type at the high (approach I) or low (approach II) resolution level. 

**2. [02.dea_MAST_statistics.R](Projects/Telomere_analysis/Single_cell_analysis/02.dea_MAST_statistics.R):** to summarize the sc-DEA with TL results at the cell type resolution level (approach I or approach II). 

**3. [Downstream_SC_DE.R](Projects/Telomere_analysis/Single_cell_analysis/Downstream_SC_DE.R):** to peform the downstream analysis on the sc-DEA with TL from the approach II.

*Of note*: The functions used in [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) and [02.dea_MAST_statistics.R](Projects/Telomere_analysis/Single_cell_analysis/02.dea_MAST_statistics.R) are defined in an [accessory script](Projects/Telomere_analysis/Single_cell_analysis/scripts/accessory_functions.R).

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
This section explains the input data and it’s structure to run the three main scripts. To follow the explanations in the **Required Input** section, you can clone the  [Projects/Telomere_analysis/Single_cell_analysis](Projects/Telomere_analysis/Single_cell_analysis) directory from the [Groningen-Microbiome repository](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) and change your current working directory.   

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

*Of note*: `git clone --filter` from git 2.19 is used to clone only objects in the [Projects/Telomere_analysis/Single_cell_analysis](Projects/Telomere_analysis/Single_cell_analysis) directory from the [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository.

If this option is not working for you, you could clone the whole [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository and change your current working directory. However, clonning the whole repository could take a larger amount of time/memory.   

```
git clone https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome
cd Groningen-Microbiome/Projects/Telomere_analysis/Single_cell_analysis/
```

### Test Data
We have provided a **cell-type-specific example dataset** for each of the **two approaches**: CD8T memory cells for approach I (n =  11,071) and CD8T cells for approach II (n =  13,246), together with some other inputs required to run the sc-DEA with TL [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R). These files are hosted in the [**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/) in the [MOLGENIS cloud](https://www.molgenis.org/).

Here, you can find the structure of the **input directory (path/to/in_dir)** for [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R). You can find the same structure in the input directory example in [**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/):

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

*Of note*: Some of the files are **strictly required**, but some of them are **optional**:

**Required**: 
* You should have a directory for each of the cell type classification levels (e.g., approach_I and approach_II directories) containing the seurat object for a specific cell type (e.g., approach_I > CD8T_memory.rds, and approach_II > CD8Tcells.rds).
* You should have a tsv separated file with the covariates considered in the scDEA model (e.g., covariates.approach_I.tab and covariates.approach_II.tab). The details of this file will be explained in the next section.

**Optional**: 
* Since the scDEA with TL is peformed at the gene-wise level across all the expressed genes, it can take a lot of time to run depending on the number of genes anc cells in your cell-type-specific seurat object. Thus, we provide several tsv files to test [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) only with a subset of genes (e.g., top3genes.approach_I.tab, top3genes.approach_II.tab, top10genes.approach_I.tab, and top10genes.approach_II.tab). The details of this file will be explained in the next section.

### Required Data 
**input directory (path/to/in_dir)**    
├── approach_I  
│   └── CD8T_memory.rds  
├── approach_II  
│   └── CD8Tcells.rds  
├── covariates.approach_II.tab  
├── covariates.approach_I.tab    

#### Cell-type-specific seurat object
A subset of the whole sc-RNAseq dataset which corresponds to a **seurat object** for a **specific cell type**. It should be located in a directory for each of the cell type classification levels (e.g., *approach_I* and *approach_II* directories). 

*Of note*:
* It must be a seurat object saved in .rds format, and named as ${cell_type}.rds (e.g., *approach_I > CD8T_memory.rds*, and *approach_II > CD8Tcells.rds*).

#### Covariates information file
A tsv file that has in the:
* 1st column: Covariates name. 
* 2nd column: Covariates type. 

*Of note*:
* Tsv separated.
* This file must have this header. 
* It is assumed that the covariates names are columns of the metadata slot of the seurat object (except from *cngeneson* which is internally calculated in the [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) script). 
* The covariates information files provided for the test datasets are the following:

* Approach I: *covariates.approach_I.tab*

| covariate  | type | 
| ------------- | ------------- | 
| TL  | contrast  |
| cngeneson | fixed  | 
| Sex | fixed  | 
| Age | fixed  | 
| donor | random  | 
| lane | random  | 

* Approach II: *covariates.approach_II.tab*

| covariate  | type | 
| ------------- | ------------- | 
| TL  | contrast  |
| cngeneson | fixed  | 
| Sex | fixed  | 
| Age | fixed  | 
| celltype.l2 | fixed  | 
| donor | random  | 
| lane | random  | 

*Of note*:
* Biological variables:
*TL* = telomere length
*Sex, Age* = donor metadata
*donor* = donor ID
*celltype.l2* = high level cell-type classification (e.g., Azimuth's prediction). You can find a detailed explanation in the methods of the [manuscript](https://www.biorxiv.org/content/10.1101/2021.12.14.472541v1).

* Technical variables:
*cngeneson* = cellular detection rate (recommended by the [MAST authors](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5))
*lane* = sequencing lane (batch effect)

### Optional Data
**input directory (path/to/in_dir)**    
├── top10genes.approach_II.tab  
├── top10genes.approach_I.tab  
├── top3genes.approach_II.tab  
└── top3genes.approach_I.tab

#### Top N single-cell differentially expressed genes (DEGs)
Since the sc-DEA with TL is peformed at the gene-wise level across all the expressed genes, it can take a lot of time to run depending on the number of cells in your cell-type-specific seurat object.

As an example, we provide several tsv files to try [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) only with a subset of genes (e.g., *top3genes.approach_I.tab, top3genes.approach_II.tab, top10genes.approach_I.tab, and top10genes.approach_II.tab*). In this case, we are picking the top 3 or top 10 DEGs in CD8T memory cells (approach I) or CD8T cells (approach II) found in the supplementary tables in the manuscript ([Table_S7.1 and Table S_7.2](https://www.biorxiv.org/content/biorxiv/early/2021/12/15/2021.12.14.472541/DC1/embed/media-1.xlsx?download=true)).

A tsv file that has in the:
* 1st column: Gene names. 

*Of note*:
* Tsv separated.
* This file do not have to incorporte a header.
* It is assumed that the gene names are in the seurat object after the pre-processing gene-level filtering performed in [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R)
* The top N DEGs files provided for the test datasets are the following:

* Approach I: Top 3 (*top3genes.approach_I.tab*) or top 10 (*top10genes.approach_I.tab*) DEGs in CD8T memory cells. See the supplementary table S7.1 in the manuscript ([Table_S7.1](https://www.biorxiv.org/content/biorxiv/early/2021/12/15/2021.12.14.472541/DC1/embed/media-1.xlsx?download=true)). For example, *top3genes.approach_I.tab*:

|   |  | 
| ------------- | ------------- | 
| DNAJA1  |  
| MTFP1 |  
| CDC42SE1 |   

* Approach II: Top 3 (*top3genes.approach_II.tab*) or top 10 (*top10genes.approach_II.tab*) DEGs in CD8T memory cells. See the supplementary table S7.2 in the manuscript ([Table_S7.2](https://www.biorxiv.org/content/biorxiv/early/2021/12/15/2021.12.14.472541/DC1/embed/media-1.xlsx?download=true)). For example, (*top3genes.approach_II.tab*):

|   |  | 
| ------------- | ------------- | 
| TMSB10  |  
| GNLY |  
| TAGAP |   


## Running the sc-DEA with TL
*Of note*: The functions used in [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) and [02.dea_MAST_statistics.R](Projects/Telomere_analysis/Single_cell_analysis/02.dea_MAST_statistics.R)are defined in an [accessory script](Projects/Telomere_analysis/Single_cell_analysis/scripts/accessory_functions.R).

**1. Clone this Github repository:** If you have not done it yet, the first step would be to clone the objects in the [Projects/Telomere_analysis/Single_cell_analysis](Projects/Telomere_analysis/Single_cell_analysis) directory from the [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository.

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

*Of note*: `git clone --filter` from git 2.19 is used to clone only objects in the [Projects/Telomere_analysis/Single_cell_analysis](Projects/Telomere_analysis/Single_cell_analysis) directory from the [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository.

If this option is not working for you, you could clone the whole [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository and change your current working directory. However, clonning the whole repository could take a larger amount of time/memory.   

```
git clone https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome
cd Groningen-Microbiome/Projects/Telomere_analysis/Single_cell_analysis/
```

**2. Set required environmental variables:**
```
input_directory=path/to/in_dir
output_directory=/path/to/out_dir
```

**2.1.** Approach I variables:

```
cell_level=approach_I  
cell_type=CD8T_memory
covariates_file=covariates.approach_I.tab
```

**2.2.** Approach II variables:

```
cell_level=approach_II 
cell_type=CD8Tcells
covariates_file=covariates.approach_II.tab
```

**2. Set optional environmental variables:**

Since the sc-DEA with TL is peformed at the gene-wise level across all the expressed genes, it can take a lot of time to run depending on the number of genes and cells in your cell-type-specific seurat object. As a test, you can run the sc-DEA with TL only a specific subset of genes for both of the approaches. For further details you can revisit the previous section "Top N single-cell differentially expressed genes (DEGs)".

```
genes_subset=top3genes
```

*Of note*: 
* If you want to use our testing dataset, you **do not need** to specify the `input_directory` variable. In this case, it will be the [**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/) in the [MOLGENIS cloud](https://www.molgenis.org/). 
* The subset genes filename must be: $genes_subset.$cell_level.tab (e.g., *top3genes.approach_I.tab or top3genes.approach_II.tab*)
* To run an example with our testing datasets (provided in [**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/)), please follow **section 3**.

**3. Run the [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) and [02.dea_MAST_statistics.R](Projects/Telomere_analysis/Single_cell_analysis/02.dea_MAST_statistics.R) scripts:**

As a **testing example**, we will run the sc-DEA with TL only for the top 3 DEGs for both of the approaches. You could also try to run it using the top 10 DEGs, or all the genes in the seurat object (it will take a lot of time/memory resources).

```
output_directory=out_dir.top3genes
genes_subset=top3genes
```

*Of note*:

* To directly use our testing datasets (in the cloud, downloaded when running the [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) script), you need to have **internet connection** to access data located in the MOLGENIS cloud ([**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/)). In this case, the `input_directory` will be the default one (https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/):

```
input_directory=https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/
```

* Otherwise, you could first download the data we provided in the MOLGENIS cloud ([**combio_andreu_2022** directory](https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/), and then use our testing datasets (local, already downloaded) when running the ([01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) script **without internet connection**. In this case, you will also need to define the `input_directory`:

```
wget -m -p -E -k -K -np https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/
mv downloads.molgeniscloud.org/downloads/combio_andreu_2022/* .
rm -r downloads.molgeniscloud.org
input_directory=$PWD
```

**3.1.** Perform the sc-DEA with TL per cell type at the high (approach I) or low (approach II) resolution level:

* Approach I:
```
cell_level=approach_I  
cell_type=CD8T_memory
covariates_file=covariates.approach_I.tab
Rscript 01.dea_MAST_glmer_TL.R --in_dir $input_directory --cell_level $cell_level --cell_type $cell_type --covs $covariates_file --genes $genes_subset --out_dir $output_directory
```

* Approach II:
```
cell_level=approach_II  
cell_type=CD8Tcells
covariates_file=covariates.approach_II.tab
Rscript 01.dea_MAST_glmer_TL.R --in_dir $input_directory --cell_level $cell_level --cell_type $cell_type --covs $covariates_file --genes $genes_subset --out_dir $output_directory
```

**3.2.** Summarizing the sc-DEA with TL results at the cell type resolution level (approach I or approach II):

* Approach I:
```
cell_level=approach_I  
Rscript 02.dea_MAST_statistics.R --in_dir $output_directory --cell_level $cell_level
```

* Approach II:
```
cell_level=approach_II  
Rscript 02.dea_MAST_statistics.R --in_dir $output_directory --cell_level $cell_level
```

*Of note*:

The outputs for each of the parameters settings are:  
**out_dir.top3genes/**  
├── approach_I      
│   ├── CD8T_memory    
│   │   ├── de_glmer.nagq0.rds  
│   │   ├── de_glmer_stats.nagq0.rds  
│   │   └── pbmc_sca.rds  
│   ├── dea.list_formatted.nagq0.rds  
│   └── Table_S7.1.csv  
└── approach_II  
    ├── CD8Tcells  
    │   ├── de_glmer.nagq0.rds  
    │   ├── de_glmer_stats.nagq0.rds  
    │   └── pbmc_sca.rds  
    ├── dea.list_formatted.nagq0.rds  
    └── Table_S7.2.csv  
 
* By default, the [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) script set the `--nagq0` argument to TRUE. In brief, it mitigate model fails to converge by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0)). See our previous [MAST github issue](https://github.com/RGLab/MAST/issues/148) and the [lme4 package documentation](https://www.rdocumentation.org/packages/lme4/versions/1.1-25/topics/glmer) for further details. In this case, the output filenames will contain the 'nagq0' tag. 
* By default, the [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R) script set the `--freq` argument to 0.1. It means that we will only test genes that are found in at least 0.1 of the cells since we won’t have any power to conclude much about less frequent genes (recommended by the [MAST authors](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0844-5))
* The `output_directory` will contain subidrectories for each of the cell type classification levels (e.g., approach_I and approach_II directories) containing for each cell type (e.g., approach_I > CD8T_memory or approach_II > CD8Tcells) the outputs from *3.1*: pbmc_sca.rds  and de_glmer.nagq0.rds (intermediate files) and de_glmer_stats.nagq0.rds (final output)
* The `output_directory` will contain subidrectories for each of the cell type classification levels (e.g., approach_I and approach_II directories) with the outputs from *3.2*: dea.list_formatted.nagq0.rds (used as input in [Downstream_SC_DE.R](Projects/Telomere_analysis/Single_cell_analysis/Downstream_SC_DE.R)) and Table_S7.1.csv/Table_S7.2.csv ([supplementary tables S7.1 and S7.2 in the manuscript](https://www.biorxiv.org/content/biorxiv/early/2021/12/15/2021.12.14.472541/DC1/embed/media-1.xlsx?download=true)). 


# Example outputs
We have provided the output directory for the test data in a tar.gz format:

** [out_dir.top3genes.tar.gz](Projects/Telomere_analysis/Single_cell_analysis/out_dir.top3genes.tar.gz): Outputs from running the commands in **section 3.** (*01.dea_MAST_glmer_TL.R* and *02.dea_MAST_statistics.R* scripts):
You can decompress them by:
```
tar -xvf out_dir.top3genes.tar.gz
```

## Running time and memory requirements

* [01.dea_MAST_glmer_TL.R](Projects/Telomere_analysis/Single_cell_analysis/01.dea_MAST_glmer_TL.R): To speed up the running time and improve the memory requirements when you use **all the genes** (not only a subset using the parameter `genes_subset`), we recommend to submit each of the commands in **section 3** as an independent job on your HPC infrastructure (i.e., run each job as an element of a job array). The running time and memory requirements will depend on:  
1. The **number of cells** in your dataset. In this case, CD8T memory cells for approach I (n =  11,071) and CD8T cells for approach II (n =  13,246).
2. The **number of genes** tested, after applying some internal filters by the same script. In this case, CD8T memory cells for approach I (n =  2,343) and CD8T cells for approach II (n = 2,378).

* [02.dea_MAST_statistics.R](Projects/Telomere_analysis/Single_cell_analysis/02.dea_MAST_statistics.R): The time and memory requirements to run this script are minimal. It is used to extract the files for: 1) [Downstream_SC_DE.R](Projects/Telomere_analysis/Single_cell_analysis/Downstream_SC_DE.R) (`$cell_level` > *dea.list_formatted.nagq0.rds*), and 2) [supplementary tables S7.1 and S7.2 in the manuscript](https://www.biorxiv.org/content/biorxiv/early/2021/12/15/2021.12.14.472541/DC1/embed/media-1.xlsx?download=true) (approach_I > *Table_S7.1.csv* and approach_II > *Table_S7.2.csv*).
