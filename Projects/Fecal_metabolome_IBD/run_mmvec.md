Installing and running mmvec  standalone version
===

Warning: Needs python 3.7 because it depends on tensorflow 1.14, newer versions are not compatible

```{sh}
#bash/zsh terminal

conda install python=3.7.0
conda create --name test_env python=3.7
conda activate test_env
conda install mmvec -c conda-forge 
conda install -c conda-forge tensorflow==1.14  
```

Convert txt to biom, just in case, check that all the columns are in the same order

```{sh}
biom convert -i ./taxa_IBD.txt -o taxa_IBD.biom --table-type="OTU table" --to-json   
biom convert -i ./taxa_controls.txt -o taxa_controls.biom --table-type="OTU table" --to-json
biom convert -i ./mtb_control.txt -o mtb_controls.biom --table-type="OTU table" --to-json 
biom convert -i ./mtb_IBD.txt -o mtb_IBD.biom --table-type="OTU table" --to-json     
```

Run mmvec
---

```{sh}
mmvec paired-omics --microbe-file taxa_all.biom --metabolite-file mtb_all.biom --num-testing-examples 45 --epochs 2000 --learning-rate 1e-5 --input-prior 1 --output-prior 1 --beta1 0.9 --beta2 0.95 --checkpoint-interval 12 --summary-interval 12 --batch-size 300 --summary-dir All_poch2000_batch_300
```

To check model parameters go to output log folder and check tensorflow 

```{sh}
 tensorboard --logdir .
```

Open browser and paste something like: http://arnaus-macbook-pro.localhost:6006/


Convert ranks to probabilities
---

```{py}
import pandas as pd
from skbio.stats.composition import clr_inv as softmax
ranks = pd.read_table('./IBD/latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ranks.txt', index_col=0)
probs = ranks.apply(softmax)
probs.to_csv('./IBD/IBDconditional_probs.txt', sep='\t')
```

Get ordinations that can be imported to R
---
```{sh}
conda activate qiime2-2020.11  
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
qiime tools import --input-path ./latent_dim_3_input_prior_1.00_output_prior_1.00_beta1_0.90_beta2_0.95_ordination.txt --output-path ./biplot.qza --type "PCoAResults % Properties('biplot')"
```


Installing and running mmvec QIIME version version
===


Install
---

```{sh}
conda activate qiime2-2020.11  
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
conda install mmvec -c conda-forge
```

Run
---
```{sh}
qiime tools import --input-path taxa_all.biom --output-path taxa_all.qza --type "FeatureTable[Frequency]" --input-format BIOMV100Format
qiime tools import  --input-path mtb_all.biom  --output-path mtb_all.qza --type "FeatureTable[Frequency]" --input-format BIOMV100Format

qiime mmvec paired-omics --i-microbes taxa_all.qza --i-metabolites mtb_all.qza --p-num-testing-examples 45 --p-epochs 2000 --p-learning-rate 1e-5 --p-input-prior 1 --p-output-prior 1  --p-summary-interval 12 --p-batch-size 300 --o-conditionals ./conditionals.qza --o-conditional-biplot ./biplot.qza --output-dir mmvec_in_qiime 
```

Visualization via Emperor
----

```{sh}
qiime metadata tabulate --m-input-file ./conditionals.qza --o-visualization conditionals-viz.qzv
qiime emperor biplot --i-biplot biplot.qza --o-visualization emperor_n100.qzv --m-sample-metadata-file metadata_mtb.txt --p-number-of-features 100
```


Plot using R
===

```{r}
library(tidyverse)
library(qiime2R)
library (ggplot2)
library(viridis)
library (plotly)
library (reshape2)

#Import ordenation file
all_ord=read_qza("./biplot.qza")

#Import metabolite info
info <- read.delim("/metadata_mtb.txt")

colnames(info)[1]="SampleID"


baseplot_cnt = ggplot() + 
  geom_point(
    data=all_ord$data$Vectors%>%
      left_join(info),
    aes(x=PC1, y=PC2, fill=SIG), alpha=0.85, 
    shape=21, size = 2.5
  ) + scale_fill_manual(values=c("white", "orange1"))

baseplot_cnt +
    theme_bw() +
    xlab(paste(round(100*all_ord$data$ProportionExplained[1],2),"%")) +
    ylab(paste(round(100*all_ord$data$ProportionExplained[2],2),"%")) +
    ggtitle("Biplot mmvec") +
    geom_segment(size=2, data=all_ord$data$Species %>% 
                     mutate(a=sqrt(PC1^2+PC2^2)) %>% 
                     top_n(10, a) %>% #keep 8 furthest away points
                     mutate(PC1=PC1*0.003, PC2=PC2*0.003), 
                 aes(x=0, xend=PC1, y=0, yend=PC2,color=FeatureID),
                 arrow = arrow(length = unit(0.3,"cm"))
    ) + scale_color_manual(values =c ("gray70","red3" ,"black","red","magenta", "steelblue2","gray40" ,"lightblue","gray80", "violetred" ,"purple"))

+ theme (legend.position = "none")+ scale_color_manual(values =c ("pink3", "salmon3", "red3", "violetred", "firebrick2", "pink", "tomato1","lightblue2", "salmon", "red" ))

```
