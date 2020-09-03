from pathlib import Path

#LL score
dic_ID = {}
with open("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/lifelines_diet_score/Lifelines_Diet_Score_InclQuintile.dat") as F:
	for line in F:
		if "PSEUDOIDEXT" in line: continue
		l = line.rstrip().split()
		LL_diet = l[1]
		dic_ID[l[0]] = {"LL" : l[1]}
print(len(dic_ID.keys()))
Column_name = ["Participant", "LL"]
for File in Path("TABLES/").glob("*_score.tsv"):
	Column_name.append(File.name.split("_")[0])
	print(File.name.split("_")[0])
	with open(File) as f:
		for line in f:
			if "Participant" in line: continue
			l =  line.rstrip().split()
			if l[0] not in dic_ID: dic_ID[l[0]] = {}
			dic_ID[l[0]][File.name.split("_")[0]] = l[1]


with open("TABLES/DIET_SCORES.tsv","w") as O:
	O.write("\t".join(Column_name) + "\n")
	for ID in dic_ID:
		INPUT = ID
		for Key in Column_name[1:]:
			if Key in dic_ID[ID]: ADD = dic_ID[ID][Key]
			else: ADD= "NA"
			INPUT += "\t" + ADD
		INPUT += "\n"
		
		O.write(INPUT)
