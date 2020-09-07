ml RPlus
echo "Running Mediterranean score script"
Rscript Create_scores/aMED_Score_LL_Cluster.R
echo "Running HEI score script"
Rscript Create_scores/HEI_Score_LL_Cluster.R
echo "Running Protein score script"
Rscript Create_scores/Protein_Score_LL_Cluster.R
echo "Running mMED script (based on Kalilli et al.)"
Rscript Create_scores/mMED_Khalilli.R

