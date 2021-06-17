
library(tidyverse)
myPaths <- .libPaths()
.libPaths(c("packages/", myPaths[1], myPaths[2]))
library(TwoSampleMR)

args = commandArgs(trailingOnly=TRUE)
args[1] -> Exposure
args[2] -> Outcome
args[3] -> Dir

print("Reading exposure as tibble")
Info_exp = read.table(Exposure, header=T) %>% as_tibble()
print("Reading exposure as Exposure data")
exposure_dat  = read_exposure_data(filename = Exposure,
    sep = "\t",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "eff.allele",
    other_allele_col = "ref.allele",
    #eaf_col = "a1_freq",
    pval_col = "P.weightedSumZ",
    #units_col = "Units",
    #gene_col = "Gene",
    samplesize_col = "N"
)
print("Getting subset of chromosomes to check in outcome")
unique(exposure_dat$chr.exposure) -> chr_interest
print("Clumping exposure data")
exposure_dat <- clump_data(exposure_dat, 
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR")

print("Chaging exposure rsID by chr:pos so that it matches outcome")
Info_exp %>% mutate(SNP = rsID) %>% select(SNP, `chr.bp`) ->Info_exp
left_join(exposure_dat, Info_exp, by="SNP") %>% mutate(SNP=`chr.bp`) %>% select(-`chr.bp`) -> exposure_dat
print(as_tibble(exposure_dat))

#MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Effect  StdErr  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal Ntotal
print("Reading outcome, this one takes some time")
outcome_dat <- read_outcome_data(Outcome,
    sep = "\t",
    snp_col = "MarkerName",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value",
    #units_col = "Units",
    #gene_col = "Gene",
    samplesize_col = "Ntotal"
)
print("Removing SNP nucleotides from SNP name and keeping chromosomes of interest")
Variants  = outcome_dat$SNP
lapply(Variants, FUN= function(x){ str_split(x, ":")[[1]][1] }) %>% as_vector() -> CHR
outcome_dat %>% mutate(CHR = CHR) %>% filter(CHR %in% chr_interest) %>% select(-CHR) -> outcome_dat
Variants  = outcome_dat$SNP
lapply(Variants, FUN= function(x){ str_split(x, "_")[[1]][1] }) %>% as_vector() -> Variants
outcome_dat$SNP = Variants

#print(Variants)
#Variants_out = c()
#for (V in Variants){
# print(V)
# Variants_out = c(Variants_out, ieugwasr::variants_to_rsid(V))
# print( ieugwasr::variants_to_rsid(V))
#}
#Variants_out = ieugwasr::variants_to_rsid(Variants)
#outcome_dat$SNP = Variants_out

print("Data harmonisation")
dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat
)

######################################
########Performing tests##############
######################################

Defaults = mr_method_list() %>% as_tibble() %>% filter(use_by_default==T) 
print("Running mendelian randomization") #Default tests, see mr_method_list() %>% as_tibble() %>% filter(use_by_default==T) 
res <- mr(dat, method_list = Defaults$obj)
print(res)

#Heterogeneity test
print("Running heterogeneity test")
Hetero_test = mr_heterogeneity(dat)
print(Hetero_test)

#Pleiotrophy
print("Pleiotropy check using mrEgger (intercept != 0?) ")
Pleiotropy_check = mr_pleiotropy_test(dat)
print(Pleiotropy_check)

#SNP individual
Sin<- mr_singlesnp(dat, all_method=c(Defaults$obj))

#Plots
print("Plotting results")
p1 <- mr_scatter_plot(res, dat) + theme_bw()
p2 <- mr_forest_plot(Sin)



print("Saving results")
paste(c(Dir, "MR_results.tsv"), collapse="") -> MR_results
paste(c(Dir, "HEterogeneity_check.tsv"), collapse="") -> Het_check
paste(c(Dir, "Pleiotropy_check.tsv"), collapse="") -> Pleiotropy_check_f
paste(c(Dir, "Individual_SNP.tsv"), collapse="") -> SNP_check_f
paste(c(Dir, "MR_plot.pdf"), collapse="") -> Plot
paste(c(Dir, "SNP_plot.pdf"), collapse="") -> Plot2

print("writing tables")
write_tsv(res, MR_results) ; write_tsv(Hetero_test, Het_check) ;  write_tsv(Pleiotropy_check, Pleiotropy_check_f) ; write_tsv(Sin, SNP_check_f)
print("saving plots")
ggsave(Plot, p1[[1]]) ; ggsave(Plot2, p2[[1]])

print("combining results")
all_res<-combine_all_mrresults(res,Hetero_test, Pleiotropy_check, Sin, ao_slc=F, Exp=T, split.exposure=T,split.outcome=F) %>% as_tibble()
paste(c(Dir, "Summary.tsv"), collapse="") -> Summ
write_tsv(all_res, Summ)
print("HTML report")
mr_report( dat, output_path = Dir)








