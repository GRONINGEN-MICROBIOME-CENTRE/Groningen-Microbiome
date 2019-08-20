
#Creator: Arnau Vich
#Year: 2018
# This code is an implementation of METAL in R. For more information please check: https://genome.sph.umich.edu/wiki/METAL_Documentation

#Approach 1: For data with P-values and sample size do: 
#-----------------------------
# 

# P=c(0.06206254,0.4491178,0.1313594, 0.02344002)
# W=c(126,208,231,904)

## Calculate P-value for your meta-analysis (can be done with sumz function from metap package which will give you 1-sided p-value)

#Square-root of sample size
Wi=sqrt(W)
#Convert p-values to Z-scores
Zi=qnorm(1-(P/2))
# Meta-zscore  
Z=(sum(Zi*Wi)/sqrt(sum(Wi^2)))
# Convert Z-score to p-value
MetaP=2*pnorm(-abs(Z))

#Cochran Q-test (test heterogeneity of your meta-analysis)

WeightedZ= sum(sqrt(W)*Zi)
totalSample=sum(W)

#Calculate expected Z
expected_Z= sqrt(W)*WeightedZ/totalSample
het_Sum= sum((Zi-expected_Z)*(Zi-expected_Z))

#Get p-value of heterogeneity test!
my_pvalue_co=pchisq(het_Sum, lower.tail = F, df=3)




#Approach 2: For data with effect size and standard errors
#-------------------------------------------------------------

# TE = vector of effect size 
# SE = vector of standard errors 
# df = degrees of freedom -> n-1 

#calculate inverse variance per cohort
selection$inverse_var.c1=1/selection$SE.c1^2
selection$inverse_var.c2=1/selection$SE.c2^2
selection$inverse_var.c3=1/selection$SE.c3^2

#Calculate SE
selection$se=sqrt(1/(selection$inverse_var.c1+selection$inverse_var.c2+selection$inverse_var.c3))

#Calculate Beta
selection$beta=(selection$inverse_var.c1*selection$Coef.c1+selection$inverse_var.c2*selection$Coef.c2+selection$inverse_var.c3*selection$Coef.c3)/(selection$inverse_var.c1+selection$inverse_var.c2+selection$inverse_var.c3)
	
#Calculate Z-score
selection$Z=selection$beta/selection$se
	
#Calculate meta p-value
selection$P=2*pnorm(-abs(selection$Z))
	
#Adjust pvalue with FDR
selection$FDR=p.adjust(selection$P,method = "fdr")

library (meta)
het=metagen(TE,SE)
my_pvalue_co=pchisq(het$Q,df=2,lower.tail=F)
