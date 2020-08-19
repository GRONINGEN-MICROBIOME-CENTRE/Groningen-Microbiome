# Calculation IBD genetic risk score


*step.1 Extract IBD GWAS summary statistics from public database*
---

Database is from the largest IBD GWAS so far https://www.ebi.ac.uk/gwas/ (reference:https://www.nature.com/articles/ng.3760#Sec12)

ebi-a-GCST004132 is CD

ebi-a-GCST004133 is UC

ebi-a-GCST004131 is IBD

```
library(TwoSampleMR)
library(MRInstruments)

gwas.hits=read.table("PRS/snps.Wehave.txt",sep = "\t",header = T,stringsAsFactors = F)

IBD.data <- extract_outcome_data(
  snps = gwas.hits$OurVariants[gwas.hits$Disease=="IBD" | gwas.hits$Disease=="IBD,UC" | gwas.hits$Disease=="UC,IBD"],
  outcomes = 'ebi-a-GCST004131'
)
CD.data <- extract_outcome_data(
  snps = gwas.hits$OurVariants[gwas.hits$Disease=="CD" | gwas.hits$Disease=="CD,UC"],
  outcomes = 'ebi-a-GCST004132'
)
UC.data <- extract_outcome_data(
  snps = gwas.hits$OurVariants[gwas.hits$Disease=="UC" | gwas.hits$Disease=="CD,UC" | gwas.hits$Disease=="IBD,UC" | gwas.hits$Disease=="UC,IBD"],
  outcomes = 'ebi-a-GCST004133'
)

IBD.data=IBD.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome")]
CD.data=CD.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome")]
UC.data=UC.data[,c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","pval.outcome","outcome")]

write.table(IBD.data,file = "IBD.summary.txt",sep = "\t",quote = F,row.names = F)
write.table(UC.data,file = "UC.summary.txt",sep = "\t",quote = F,row.names = F)
write.table(CD.data,file = "CD.summary.txt",sep = "\t",quote = F,row.names = F)
```

*step.2 Check GWAS SNPs/SNP-proxy in our genetics dataset*
---

241 indepedent-variants in total from the largest GWAS (Europe+Asian) so far. However, 13 variants can't find back in our genetics data (WES+GSA), of which six have proxies. The rest of seven varaints or their proxy are not found due to the low allele frequency in European populations (high in Asian).

Variants are not found:rs2636006, rs111781203, rs56691573, rs75900472, rs11235667, rs7329174 and rs80244186

Variants are replaced with proxy based on LD (r=1): rs3180018, rs144344067, rs2593855, rs71559680, rs149169037 and rs200349593

The final variants number to be used in PRS calculation is 234.


*step.3 calculate PRS based on GWAS summary in our datasets*
---

Formular of PRS: sum(beta * effect allele) / total number of effect alleles

reference(https://www.prsice.info/step_by_step/; Choi,et,al.(2020).NatureProtocols)

```
Rscript ~/PRSice.R  --dir . --prsice ~/PRSice_linux --base CD.summary.txt --target /groups/umcg-gastrocol/tmp04/Metabolic_Project/Genetics/Combined/All.merged.1_22chr.biallic --stat beta.outcome --binary-target F --out CD.PRS --no-regress T --no-full T --snp SNP --chr chr --bp pos --A1 effect_allele.outcome --A2 other_allele.outcome --pvalue pval.outcome

script ~/PRSice.R  --dir . --prsice ~/PRSice_linux --base UC.summary.txt --target /groups/umcg-gastrocol/tmp04/Metabolic_Project/Genetics/Combined/All.merged.1_22chr.biallic --stat beta.outcome --binary-target F --out UC.PRS --no-regress T --no-full T --snp SNP --chr chr --bp pos --A1 effect_allele.outcome --A2 other_allele.outcome --pvalue pval.outcome

Rscript ~/PRSice.R  --dir . --prsice ~/PRSice_linux --base IBD.summary.txt --target /groups/umcg-gastrocol/tmp04/Metabolic_Project/Genetics/Combined/All.merged.1_22chr.biallic --stat beta.outcome --binary-target F --out IBD.PRS --no-regress T --no-full T --snp SNP --chr chr --bp pos --A1 effect_allele.outcome --A2 other_allele.outcome --pvalue pval.outcome
```
