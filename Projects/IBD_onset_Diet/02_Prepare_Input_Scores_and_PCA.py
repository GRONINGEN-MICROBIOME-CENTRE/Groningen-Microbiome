

dic_grouping = {}

#Groups_for_PCA  Item_name_cluster       (Groups_for_scores)     Item_name_cluster       Item_name_Eva
print("Preparing dic in which each cluster name has a category (used for PCA) and a name in Eva's scripts")
with open("TABLES/Table_groups.tsv") as F:
	for line in F:
		if "Groups_for_PCA" in line: continue
		l = line.rstrip().split("\t")
		dic_grouping[l[1]] = { "PCA" : l[0], "Diet_score" : l[4] }
		dic_grouping[l[3]] = { "PCA" : l[0], "Diet_score" : l[4] }

def Make_PCA_input(File):
	'Generate the input for PCA'
	dic_categories = {}
	dic_samples = {}
	list_samples = []
	with open(File) as F:
		n = 0
		for line in F:
			l = line.rstrip().split("\t")
			#Go food item by food item and transform
			if n == 0:
				Header = l[0:27]
				for m in range(27, len(l)):
					item  = l[m]
					if item in dic_grouping:
						Category = dic_grouping[item]["PCA"]
						if Category not in dic_categories: dic_categories[Category] = []
						dic_categories[Category].append(m)
			else:
				if len(l) != 137: continue
				list_samples.append(l[0:27])
				for Category in dic_categories:
					Total = 0
					for position in dic_categories[Category]:
						Total += float(l[position])		
					if l[0] not in dic_samples: dic_samples[l[0]] = {}
					dic_samples[l[0]][Category] = Total
			n+=1

	Categories = dic_categories.keys()
	Header.extend(Categories)
	with open(File + "_PCA_input", "w") as F:
		F.write("\t".join(Header)+"\n")
	for Sample in list_samples:
		ID = Sample[0]
		Input = Sample
		for Category in Categories:
			Count = dic_samples[ID][Category]
			Input.append(str(Count))
		with open(File + "_PCA_input", "a") as F:
			F.write("\t".join(Input)+"\n")

def Change_names_Eva_script(File):
	'Generates the input for score scripts'
	with open(File) as F:
		n=0
		Total = ""
		to_keep = []
		for line in F:
			l = line.rstrip().split("\t")
			if n == 0:
				New_header = l[0:27]
				for m in range(27, len(l)):
					item = l[m]
					if item not in dic_grouping: continue
					to_keep.append(m)
					New_header.append(dic_grouping[item]["Diet_score"])
				n +=1
				continue
			New_line = "\t".join(l[0:27])
			if len(l) != 137: continue
			for keep in to_keep:
				New_line += "\t" +  l[keep]
			Total += New_line + "\n"
	with open(File + "_script_eva_input", "w") as F:
		F.write("\t".join(New_header) + "\n")
		F.write(Total) 


for i in ["TABLES/Table_participant_GDAG.tsv","TABLES/Table_participant_KCAL.tsv"]:
	print(i)
	if i != "TABLES/Table_participant_KCAL.tsv":
		print("Generating PCA input")
		Make_PCA_input(i)
	print("Generating score script input")
	Change_names_Eva_script(i)



