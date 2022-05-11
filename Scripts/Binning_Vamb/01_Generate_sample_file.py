from pathlib import Path

#This iteration will be different depending on the data structure of the directory 

with open("samples2data.txt", "w") as O: pass
for F in Path("/groups/umcg-dag3/tmp01/HIV_microbiome/rawdata/"). glob("*_1.fq.gz"):
	N = F.name.split("_1")[0]
	print(F)
	R = str(F).replace("_1.fq.gz","_2.fq.gz")
	with open("samples2data.txt", "a") as O:
		O.write( "\t".join([N, str(F), R]) + "\n")
	
