# A simple calculation of IBD poly genetic risk score


*step.1 Extracted IBD GWAS summary statistics from public database*
---

Database is https://www.ebi.ac.uk/gwas/ (reference:https://www.nature.com/articles/ng.3760#Sec12)

ebi-a-GCST004132 is CD

ebi-a-GCST004133 is UC

ebi-a-GCST004131 is IBD


*step.2 Check GWAS SNPs/SNP-proxy in our genetics dataset*
---

241 indepedent-variants in total from the largest GWAS (Europe+Asian) so far. However, 13 variants can't find back in our genetics data (WES+GSA), of which six have proxies. The rest of seven varaints or their proxy are not found due to the low allele frequency in European populations.

Variants are not found:rs2636006, rs111781203, rs56691573, rs75900472, rs11235667, rs7329174 and rs80244186

Variants are replaced with proxy based on LD (r=1): rs3180018, rs144344067, rs2593855, rs71559680, rs149169037 and rs200349593
