

#1. Make reference set (all observed variations, add a threshold of number of samples it needs to be observed)
#2. Given the calls per sample, add 1 to the positiosn shared with reference and 0 to the positions not shared
#3. Compute genetic distance between profiles from 2. Manhattan distance and Fst --> https://github.com/sergioSEa/lonely_scripts/blob/master/Renelies-Hamilton2019/SNV_filtering.py

from pathlib import Path
import sys
import numpy as np
from scipy import spatial,stats,cluster
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
print("> Running script")


def Find_variants(Tool):
	dic_overall_variation = {}
	Sample_list = []
	Pattern = "SRS*/Calls/{Tool}.tsv".format(Tool=Tool)
	for F in Path(".").glob(Pattern):
		Sample = str(F).split("/")[0]
		Sample_list.append(Sample)
		with open(F) as Sample_dir:
			for line in Sample_dir:
				l = line.rstrip().split()
				ID = "-".join([l[0], l[1], l[3]])
				if ID not in dic_overall_variation:
					dic_overall_variation[ID]= []
				dic_overall_variation[ID].append(Sample)
	return(dic_overall_variation, Sample_list)
def Create_variation_profile(Tool,Threshold, Tool2="No"):
	print("> Starting to compile variants")
	dic_overall_variation = {}
	Sample_list = []
	
	dic_overall_variation, Sample_list = Find_variants(Tool)
	if Tool2 != "No":
		dic_overall_variation2, Sample_list = Find_variants(Tool2)
		dic_overall_variation_common= {}
		for Variant in dic_overall_variation2:
			if Variant not in dic_overall_variation: continue
			Samples_2 = dic_overall_variation2[Variant] ; Samples_1 = dic_overall_variation[Variant]
			Common = list(set(Samples_1) & set(Samples_2))
			dic_overall_variation_common[Variant] = Common
		dic_overall_variation = dic_overall_variation_common
	#If filter set, remove overall variation with fewer samples than threshold
	if Threshold != 0:
		print("Applying thresholds on variant population prevalence")
		Reference_variation = {}
		for Variant in dic_overall_variation:
			N = len(dic_overall_variation[Variant])
			if N < Threshold: continue
			Reference_variation[Variant] = dic_overall_variation[Variant]

	else: Reference_variation = dic_overall_variation
	print("Recovered {x} variants from {y} samples ; under prevalence threshold {z} variants remained".format(x=len(dic_overall_variation.keys()), y = len(Sample_list), z=len(Reference_variation.keys())))
	return(Reference_variation, Sample_list)
def Create_variation_profile_taxa(Tool,Variants,Taxa, Threshold=0):
	print("> Filtering variants for  taxa")
	Reference_variation = {}
	for Variant in Variants:
		CHR = Variant.split("-")[0]
		if CHR in Taxa:
			N = len(Variants[Variant])
			if N < Threshold: continue
			Reference_variation[Variant] = Variants[Variant]
	return(Reference_variation)

				
def SFS_plot(Tool, Plot_name ,Taxa= "No"):
	print("> Starting to compile variants for SFS")
	dic_overall_variation, Sample_list = Find_variants(Tool)
	print("> Computing frequency per allele")
	Variation_freq = {}
	Total_number_variants = len(dic_overall_variation.keys())
	for Variant in dic_overall_variation:
		if Taxa !=  "No": 
			CHR = Variant.split("-")[0]
			if CHR not in Taxa: continue
		N = len(dic_overall_variation[Variant])
		if N not in Variation_freq: Variation_freq[N] = 0
		Variation_freq[N] += 1
	data = list(Variation_freq.items())
	data = np.array(data)
	print("Plotting  site frequency spectrum (SFS)")
	
	rects1 = plt.bar(data[:,0], data[:,1])
	print("Total number of variants = {N}".format(N=str(Total_number_variants)))
	plt.xlabel('Derived allele frequency')
	plt.ylabel('Number of sites')
	labels = data[:,0]
	plt.xlim([1, 86])
	plt.savefig(Plot_name,bbox_inches='tight')
	plt.clf() 	
def Sample_variants(Reference_variation, Sample_list):
	print("> Creating a binary profile of reference variants per sample")
	#Create profile per sample
	Reference_size = len(Reference_variation.keys())
	Sample_profile = {}
	for S in Sample_list:
		Sample_profile[S] =  [0] * Reference_size

	Position = 0
	for Variant in Reference_variation:
		for Sample  in Reference_variation[Variant]:
			Sample_profile[Sample][Position] = 1
		Position += 1
	print("Sample profile is finished")
	
	return(Sample_profile)


def Distance(Sample_profile):
	print("> Computing genetic distance between samples based on each sample's SNV profile")
	#Compute distnace between samples
	Matrix_distance = np.zeros((len(Sample_profile.keys()),len(Sample_profile.keys())))
	Input_matrix =  False
	index1= 0
	print("Preparing matrix")
	for item in Sample_profile:
		Profile_1 = np.array(Sample_profile[item])
		Profile_1 = Profile_1.reshape(1,Profile_1.shape[0])
		if isinstance(Input_matrix, (bool)): Input_matrix= Profile_1
		else: Input_matrix = np.vstack((Input_matrix, Profile_1))
	print(Input_matrix.shape)
	print("Computing manhattan distance")
	Matrix_distance = spatial.distance.pdist(Input_matrix, metric='cityblock')	
	#print(Matrix_distance.shape)
	print("Distance matrix computed, shape:")
	print(spatial.distance.squareform(Matrix_distance).shape)	
	return(Matrix_distance)

def Dendogram_plotting(Cluster,Sample_list, Index_list,Dedogram_name):
	print("Generating dendogram")
	plt.figure()
	dn = cluster.hierarchy.dendrogram(Cluster,labels=Sample_list)
	#label_colors = {'a': 'r', 'b': 'g', 'c': 'b', 'd': 'm'}
	#ax = plt.gca()
	#xlbls = ax.get_xmajorticklabels()
	#for lbl in xlbls:
	#	lbl.set_color(label_colors[lbl.get_text()])
	plt.savefig(Dedogram_name)
	plt.clf()

def Clustering_samples(Distance_matrix,Sample_list, Index_list,Dedogram_name):
	print("> Using distance matrix for clustering")
	Cluster = cluster.hierarchy.linkage(y=Distance_matrix, method="single") #Nearest point algorithm, distance between clusters in the minimal distance between its members | y needs to be a condensed distance matrix, as output from scipy dist
	print("Cluster completed - Linkage matrix represents 1.Index sample, 2.Index sample, 3.distance, 4.number samples in file; each row is a cluster, once cluster is done, the index of the cluster is max(Index_sample)+row_cluster ")
	Dendogram_plotting(Cluster,Sample_list, Index_list,Dedogram_name)
	print("> Count how many clusters are between baseline and followup from the same sample")
	Identical = 0
	Not_identical = 0
	for HMP_sample in Index_list:
		#For eachHMP sample we have the index, thus, we take the row where this index is found. We first try in the first column, if no match,, try in the second 
		Indexes  = Index_list[HMP_sample]
		mask = Cluster[:,0] == Indexes[0]
		Row = Cluster[mask, :]
		if Row.shape[0] == 0: Row = Cluster[Cluster[:,1]==Indexes[0], :]
		if Row[0,0] in Indexes  and  Row[0,1] in  Indexes: Identical+=1
		else: Not_identical+=1
	print("Number of participants were samples in basiline and follow up cluster together: {I}\nNumber of participants that do not: {N}".format(I=str(Identical), N=str(Not_identical)))
	return(str(float(Identical)/(Identical+Not_identical)))

def Make_index_list(Sample_list, HMP):
	#Make a list of indexes from Sample_list that go together in HMP
	Index_list = {}
	Index = 0
	for item in Sample_list:
		if HMP[item] not in Index_list:
			Index_list[HMP[item]] = []
		Index_list[HMP[item]].append(float(Index))
		Index += 1
	return(Index_list)

def Compare_distances(Distance_matrix,Index_list, Boxplot_name):
	print("> Comparing genetic distances beween same individual and different samples")
	Distance_matrix = spatial.distance.squareform(Distance_matrix)
	#Split matrix in two vectors, distances between same ID and  interindividual distance, compare by a stat test
	Intra_individual = []
	Inter_individual = []
	Intra_samples = []
	for Sample  in Index_list:
		Index = Index_list[Sample]
		Intra_samples.extend([str(int(Index[0])) + "-" +  str(int(Index[1])), str(int(Index[1])) + "-"+ str(int(Index[0]))])
		D = Distance_matrix[int(Index[0])][int(Index[1])]
		Intra_individual.append(D)

	Distance_triangular = np.triu(Distance_matrix)
	i1= 0

	for Row in Distance_triangular:
		i2 =0
		for column in Row:
			Position =  str(i1) +"-" + str(i2)
			if Position in Intra_samples: 
				i2 +=1 ; continue
			if i1 == i2:
				i2 +=1 ; continue
			if column == 0: i2 +=1 ; continue
			Inter_individual.append(column)			
			i2 +=1
		i1+=1
		
	print("Performing wilcox test between inter and intra individual distanes")
	w, p = stats.mannwhitneyu(Intra_individual, Inter_individual)
	print("Pvalue=" + str(p))
	M = abs(np.mean(np.array(Intra_individual)) - np.mean(np.array(Inter_individual)))
	print("Mean difference=" + str(M))
	print("Plotting distributions")
	data = [np.array(Intra_individual), np.array(Inter_individual)]
	fig, ax = plt.subplots()
	ax.boxplot(data)
	positions = (1, 2)
	labels = ("Intra-individual", "Inter-individual")
	plt.xticks(positions, labels)
	plt.savefig(Boxplot_name)	
	plt.clf()
def Prepare_table(Percentage, Tool, N, Clustered, Taxa="No", Thresh=1):
	if Taxa=="No":
		return(["Total SNV profile", Tool, str(N), Clustered , "NA", "NA", Thresh])
	else:
		O = [Taxa, Tool, str(N), Clustered]
		with open("info/SupTable3.tsv") as F:
			#TP      TN      FP      FN      Taxa    Tool    Length_genome   Benchmark_Design
			for line in F:
				if "TP" in line: continue
				l = line.rstrip().split()
				if l[4] != Taxa: continue
				if l[-3] != Tool: continue
				TP = float(l[0]) ; FN =float(l[3]) ; FP=float(l[2])
				SN = TP/(TP+FN)
				Precision=TP/(TP+FP)
				O.extend([str(SN), str(Precision), Thresh])
		return(O)








print("Getting info on participants")
HMP = {}
with open("Key_participants_HMP.tsv") as I_file:
	for line in I_file:
		l = line.rstrip().split()
		ID = l[1].rstrip("_F")
		HMP[l[0]]= ID
		#if ID not in HMP: HMP[ID] = []
		#HMP[ID].append(l[0])a





def Run_analysis(Tool, Threshold, SFS, tool= "", Tool2="No"):
	if tool == "": tool= Tool
	#Output files
	Dendogram_File = "Plots/Cluster_{Tool}_{T}.pdf".format(Tool=tool, T= str(Threshold))
	Boxplot_File = "Plots/Inter-vs-Intra_{Tool}_{T}.pdf".format(Tool=tool, T= str(Threshold))
	FS_Plot = "Plots/SFS_{P}.pdf".format(P=tool)
	if SFS == True:
		SFS_plot(Tool, FS_Plot)
		
	###For all bacteria at the same time
	Index_list = False
	Table = []
	
	Reference_variation, Sample_list = Create_variation_profile(Tool,Threshold, Tool2)
	Sample_profile = Sample_variants(Reference_variation, Sample_list) 
	Matrix_distance =  Distance(Sample_profile)
	Index_list =  Make_index_list(Sample_list, HMP)
	Compare_distances(Matrix_distance,Index_list,Boxplot_File)
	Clusters  = Clustering_samples(Matrix_distance,Sample_list,Index_list,Dendogram_File)

	Table.append(Prepare_table(Clusters, Tool, len(Reference_variation), Clusters, Taxa="No",Thresh=str(Threshold)))

	##Per bacteria taxa
	print("> Starting individual taxa analysis (requires PythonPlus)")
	dic_taxa = {}
	from Bio import SeqIO
	print("Compiling contigs per reference")
	with open("Bacterial_selection.tsv") as B:
		for line in B:
			l = line.split()[0]
			with open("References/{F}.fa".format(F=l)) as handle:
				Chrs = SeqIO.parse(handle, "fasta")
				dic_taxa[l] = []
				for Chr in Chrs:
					dic_taxa[l].append(str(Chr.id))
	
	
	print("Compiling all variants")
	Variants, Sample_list = Find_variants(Tool)
	for taxa in dic_taxa:
		print("> Runing for taxa: {t}".format(t=taxa))
		Dendogram_File_S =  "Plots/Cluster_{Tool}_{T}_{S}.pdf".format(Tool=tool, T= str(Threshold), S=taxa)
		Boxplot_File_S =  "Plots/Inter-vs-Intra_{Tool}_{T}_{S}.pdf".format(Tool=tool, T= str(Threshold), S=taxa)
		SFS_Plot = "Plots/SFS_{P}_{T}.pdf".format(P=Tool, T=taxa)
		if SFS == True: 
			SFS_plot(Tool, SFS_Plot ,Taxa= dic_taxa[taxa])

		Reference_variation =Create_variation_profile_taxa(Tool,Variants,dic_taxa[taxa],Threshold)
		Sample_profile = Sample_variants(Reference_variation, Sample_list)
		Matrix_distance =  Distance(Sample_profile)
	
		if Index_list == False: Index_list =  Make_index_list(Sample_list, HMP)


		Compare_distances(Matrix_distance,Index_list,Boxplot_File_S)
		Clusters  = Clustering_samples(Matrix_distance,Sample_list,Index_list,Dendogram_File_S)
		Table.append(Prepare_table(Clusters, Tool, len(Reference_variation), Clusters, Taxa=taxa, Thresh=str(Threshold)))

	Output = "Plots/Table_stats_{T}-{N}.tsv".format(T=Tool, N=Threshold)
	print("Saving results at "+Output)
	with open(Output.format(T=Tool, N=Threshold) , "w") as O:
		Header = "\t".join(["Taxa_profiled", "Tool","Number Variants", "% Clustered", "Sensitivity", "Precision", "Threshold"]) + "\n"
		O.write(Header)
		for item in Table:
			O.write("\t".join(item) + "\n")

if sys.argv[1] == "all":
	Run_analysis("haplotypecaller", 1, True)
	Run_analysis("haplotypecaller", 2, False)
	Run_analysis("mutect2", 1, True)
	Run_analysis("mutect2", 2, False)
else:
	if len(sys.argv) < 2:
		print("Please specify the tool to use: python Calculate_genetic_distance.py [haplotypecaller|mutect2|intrain]")
		exit()
	tool = sys.argv[1] ; Tool = tool
	Tool2= "No"
	if "-" in Tool:
		T = Tool.split("-")
		Tool = T[0]
		Tool2 = T[1]
	if len(sys.argv) == 3:
		Threshold = int(sys.argv[2])
		SFS = False
	elif  len(sys.argv) > 3:
		Threshold = int(sys.argv[2])
		SFS =True
	else:
		Threshold = 0
		SFS = False
	Run_analysis(Tool, Threshold, SFS,tool, Tool2) #Tool2	

