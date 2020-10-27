#from pathlib import Path
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ref', '-r', action = 'store', dest = 'ref')
parser.add_argument('--calls', '-c', action = 'store', dest = 'calls')
parser.add_argument('--covered', '-cov', action = 'store', dest = 'covered')
parser.add_argument('--output', '-out', action = 'store', dest = 'output')
parser.add_argument('--threshold', '-tresh', action = 'store', dest = 'threshold')
args = parser.parse_args()


def Chromosomal_length(File):
	dic_lengths= {}
	#Length_Chromosome = 0
	with open(File) as F:
		for line in F:
			l = line.split()
			ID = l[0]
			if ID not in dic_lengths: dic_lengths[ID] = 0
			dic_lengths[ID] += 1
			
			#Length_Chromosome += 1
	return(dic_lengths)
def Variants_to_dic(File, Threshold):
	dic_Chr_m = {}
	dic_Chr = {}
	All_chr = []
	with open(File) as F:
		for line in F:
			l = line.rstrip().split()
			if len(l) == 5 and not l[-1].isalpha():
				QUAL = float(l[-1])
				All_chr.append(l[0])
				if QUAL < float(Threshold): continue
				Chr = l[0] ; Pos = l[1] ; Ref = l[2] ; Variant = l[3]
			elif len(l) != 6: Chr = l[0] ; Pos = l[1] ; Ref = l[2] ; Variant = l[3]
			elif len(l) == 6: Chr = l[0] ; Pos = l[1] ; Ref = l[3] ; Variant = l[4]
			Variant = "{C}:{P}:{R}:{A}".format(C=Chr, P=Pos, R=Ref, A=Variant)
			if l[-1] == "Mutation":
				if Chr not in dic_Chr_m: dic_Chr_m[Chr] = []
			if Chr not in dic_Chr: dic_Chr[Chr] = []
			dic_Chr[Chr].append(Variant)
	
	for Chr in set(All_chr):
		if Chr not in dic_Chr:
			dic_Chr[Chr] = [] 
	return(dic_Chr, dic_Chr_m)




def Estimation_False_n_True_calls(Reference_variants, Called_variants,Called_variants_mutations ,Total_positions):
	Results = []
	for Chr in Called_variants:
		TP = 0 ;TN = 0 ; FP = 0 ; FN = 0
		#TT = []
		Ref_variants =  Reference_variants[Chr]
		if not Chr in Called_variants: called = []
		else : called = Called_variants[Chr]
		for Variant in Ref_variants:
			if Variant in called: TP += 1 #; TT.append(Variant)
			elif Variant not in called: FN += 1
		#for V in called:
		#	if V not in TT:
		#		print(V)
		FP += len(called) - TP
		TN = Length_Chromosome[Chr] - (TP + FP) - FN
		R = "\t".join([str(TP),str(TN),str(FP),str(FN)])
		Results.append(R)
	return(Results)

#Steps
#1. Get the total number of "calls" and "no calls", i.e all positions covered
Coverage_file = args.covered
Length_Chromosome = Chromosomal_length(Coverage_file)

#2. Save real variants into a dictionary
threshold = args.threshold
Ref_Variants_file = args.ref
dic_Chr, dic_2 = Variants_to_dic(Ref_Variants_file, threshold)

#3. Save called variants to dictironary
Called_Variants_file = args.calls

dic_Chr_called, dic_Chr_called_mutationsIntroduced = Variants_to_dic(Called_Variants_file,threshold)

#3. Check if the called variatns have been called, and count the number of TP,FP and FN. Use the total number of positions to calculate TN
if len(dic_Chr_called) != 0:
	Call_stats = Estimation_False_n_True_calls(Reference_variants = dic_Chr, Called_variants = dic_Chr_called, Called_variants_mutations=dic_Chr_called_mutationsIntroduced  , Total_positions = Length_Chromosome)

output = args.output # Working_dir+"Statistics_calling/{tool}.{ID}.tsv"
with open(output, "w") as O:
	Header = "TP,TN,FP,FN,Taxa,Tool".split(",")
	
	if "_stats" in output:
		#Working_dir + "tmp/{tool}_{ID}_"+ str(Threshold) + "_stats"
		out = output.split("/")[-1]
		Tool = out.split("_")[0]
		Taxa = "_".join(out.split(str(threshold)+"_stat")[0].split("_")[1:-1])
	else:
		out = output.split("/")[-1]
		out = out.split(".")
		Tool = out[0]
		Taxa = out[1]
	To_write = ""
	for Chr in Call_stats:
		Chr = Chr +"\t"+Taxa + "\t" + Tool + "\n"
		To_write += Chr
	O.write("\t".join(Header) + "\n" + To_write)
