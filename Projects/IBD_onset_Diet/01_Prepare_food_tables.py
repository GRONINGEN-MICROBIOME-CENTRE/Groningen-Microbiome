#Generates input for PCA script and input for diet socre generation scripts.


#Total calories consumed per patient, one row per patient
Info_participant ="tmp/Nutrient_intake_Participant_Results_NEW.dat" #"/groups/umcg-lifelines/prm03/releases/pheno_lifelines/v1/tab_separated_labels/Nutrient_intake_Participant_Results_NEW.dat"
#Food intake of patient per item
Nutrients = "tmp/Nutrient_intake_Item_Results_NEW.dat"  #"/groups/umcg-lifelines/prm03/releases/pheno_lifelines/v1/tab_separated_labels/Nutrient_intake_Item_Results_NEW.dat"

dic_participant = {}
#Header + individual item information
Head = "SumOfkcal       SumOfkJ SumOfeiwittot   SumOfeiwitplant SumOfeiwitdier  SumOfvettot     SumOfkhtot      SumOfmodis      SumOfpolys     SumOfalcohol     SumOffree_sug   SumOfadd_sug    SumOfglucose    SumOffructose   SumOflactose    SumOfmaltose    SumOfsucrose    SumOfGlyc_ind   SumOfGlyc_load  leeftijd	GESLACHT        GEWICHT LENGTE  BMI     BMR     EI_BMR".split()
#Save participant in a dic, add all nutritional information
print("Saving nutritional information per participant")
with open(Info_participant) as F:
	for line in F:
		if "PSEUDOIDEXT" in line: continue
		l = line.rstrip().split("\t")
		if l[0] not in dic_participant: dic_participant[l[0]] = []
		dic_participant[l[0]].extend(l[2:-1])


#Add to the information of the patient KCAl per item and Grams_per_day per item. Item name is in column 3, 4 Gday, 5 Kcal
columns_interest = [3,4,5] #itemnaam        gdag    kcal
Food = []
print("Adding food nutritional info")
with open(Nutrients) as F:
	for line in F:
		if "PSEUDOIDEXT" in line: continue
		l = line.rstrip().split("\t")
		if str(type(dic_participant[l[0]][-1])) !=  "<type 'dict'>": 
			dic_food = {}
			dic_participant[l[0]].append(dic_food)
		Food.append(l[3])
		dic_participant[l[0]][-1][l[3]] = [l[5],l[4]]
Food = set(Food) 
#random participant, what we want are all food items
food_items = list(Food) #dic_participant["30453317"][-1].keys()

Head.extend(food_items)

KCAL_file = "TABLES/Table_participant_KCAL.tsv"
Grams_day = "TABLES/Table_participant_GDAG.tsv"
for I in [KCAL_file,Grams_day]:
	with open(I,"w") as F:
		F.write("Participant\t" + "\t".join(Head) + "\n")

print("Saving output")
for Participant in dic_participant:
	Input1 = []
	Input2 = []
	print(Participant)
	Input_f = [Participant]
	Input_f.extend(dic_participant[Participant][0:-1])
	Input1.extend(Input_f)
	Input2.extend(Input_f)
	for Item in food_items:
		if Item not in dic_participant[Participant][-1]: Food = ["NA", "NA"]
		else: Food = dic_participant[Participant][-1][Item]
		#print(Food, len(Input1))
		Item_1 = Food[0] ; Item_2 = Food[1]
	
		Input1.append(Item_1)
		Input2.append(Item_2)
		if len(Input1) != 137: continue
	with open(KCAL_file,"a") as F:
		F.write("\t".join(Input1) + "\n")
	with open(Grams_day,"a") as F:
		F.write("\t".join(Input2) + "\n")


		
