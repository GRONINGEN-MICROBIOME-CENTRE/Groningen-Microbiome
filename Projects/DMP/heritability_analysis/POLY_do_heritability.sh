#!/bin/bash
# ===========================================
# Notes: 
# - data file (microbiome) is _dag3_clrnorm_zix.dat
# - pedigree file is _dag3_clrnorm_zix.ped
#   - these files cannot be shared here because of risk of compromising privary of study participants
#   - data are avaliable on demand, subject to approval from UMCG ethics comittee and Lifelines biobank
# - scripts assume that POLY is installed and pedtools are installed (both in system PATH)
#   - POLY link: https://people.virginia.edu/~wc9c/poly/input.html
#   - pedtools link: https://csg.sph.umich.edu/abecasis/pedstats/download/index.html

# Heritability analysis workflow:
# ==================================================================================
# A) data preparation
#   - here predstats toolkit is used to 
#     clean and trim pedigree files (by removing unlinked participants and singletons)
# ==================================================================================
# 1) Calculate pedigree statistics, rewrite/clean pedigree
pedstats -d _dag3_clrnorm_zix.dat -p _dag3_clrnorm_zix.ped --pairs --verbose | tee ___dag3_unclean_stats_pairs
pedstats -d _dag3_clrnorm_zix.dat -p _dag3_clrnorm_zix.ped --rewritePedigree | tee ___dag3_unclean_rew

# 2) Calculate statistics on rewritten pedigree, trim pedigree
#  > rename files (pedstats default output is pedstats.dat and pedstats.ped)
mv pedstats.dat _dag3_clrnorm_zix_rew.dat
mv pedstats.ped _dag3_clrnorm_zix_rew.ped
#  > calculate statistics, trim and rewrite        
pedstats -d _dag3_clrnorm_zix_rew.dat -p _dag3_clrnorm_zix_rew.ped --pairs --verbose | tee ___dag3_rewritten_stats_pairs
pedstats -d _dag3_clrnorm_zix_rew.dat -p _dag3_clrnorm_zix_rew.ped --rewritepedigree --trim | tee ___dag3_rewritten_rewtrim

# do stats on rewritten & trimmed pedigree
#  > rename files (pedstats default output is pedstats.dat and pedstats.ped)
mv pedstats.dat _dag3_clrnorm_zix_rew_trimmed.dat
mv pedstats.ped _dag3_clrnorm_zix_rew_trimmed.ped
#  > calculate statistics            
pedstats -d _dag3_clrnorm_ziz_rew_trimmed.dat -p _dag3_clrnorm_ziz_rew_trimmed.ped --pairs --verbose | tee ___dag3_rew_trimmed_stats_pairs

# B) run POLY for cleaned & trimmed pedigree files
# =====================================================
poly -d _dag3_clrnorm_ziz_rew_trimmed.dat -p _dag3_clrnorm_ziz_rew_trimmed.ped -t __model_additive.model --output ___polyrun_dag3_clrnorm.txt | tee ___polyrun_dag3_heritability.log
