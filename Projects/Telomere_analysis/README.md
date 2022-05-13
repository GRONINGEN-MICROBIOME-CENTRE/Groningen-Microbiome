# Genetic, parental and lifestyle factors influence telomere length - Analysis scripts

Here the data analysis script used in (Info publication) can be found.

Two major analyses are included:


**1. [telomers_association.R](Projects/Telomere_analysis/telomers_association.R):** to perform the telomere exploratory analyses and telomere-metadata associations. 

**2. [Single_cell_analysis](Projects/Telomere_analysis/Single_cell_analysis/):** Includes the major scripts for the single cell transcriptomic association to telomere length. 


## Contact
If you have any questions or issues, feel free to open an issue or directly email Sergio Andreu-Sánchez (s.andreu.sanchez at umcg.nl). 


## Required Software
* **R** >=4.0.0 version: You need to install the packages loaded in the [analaysis script](telomers_association.R).



## Required Input
This section explains the input data. To follow the explanations in the **Required Input** section, you can clone the  [Projects/Telomere_analysis/Single_cell_analysis](Projects/Telomere_analysis/) directory.  

```
git clone \
  --depth 2  \
  --filter=blob:none  \
  --sparse \
  https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome \
;
cd Groningen-Microbiome
git sparse-checkout set Projects/Telomere_analysis/
cd Projects/Telomere_analysis/Single_cell_analysis/
```

*Of note*: `git clone --filter` from git 2.19 is used to clone only objects in the [Projects/Telomere_analysis/](Projects/Telomere_analysis/) directory from the [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository.

If this option is not working for you, you could clone the whole [Groningen-Microbiome](https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome) repository and change your current working directory. However, clonning the whole repository could take a larger amount of time/memory.   

```
git clone https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome
cd Groningen-Microbiome/Projects/Telomere_analysis/Single_cell_analysis/
```

**Required**: 
* To run the script you will need to get access to phenotypic information for Lifelines participants. Please contact Lifelines if interested.