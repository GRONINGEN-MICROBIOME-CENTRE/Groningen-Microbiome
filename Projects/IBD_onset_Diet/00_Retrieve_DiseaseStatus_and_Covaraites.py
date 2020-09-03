dic_covariates = {}

Path = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/lifelines_diet_score/"
Possible_files = ["Participant (Pa_99_G).dat", "Measurement_Anthropometry.dat", "Questionnaire_Smoking.dat", "Nutrient_intake_Participant_Results.dat"]

for File in Possible_files:
        #File = "Nutrient_intake_Participant_Results.dat"
        To_open = "/groups/umcg-lifelines/prm03/releases/pheno_lifelines/v1/tab_separated_labels/"+File
        dic_entries = {"Age_n":False, "Sex_n": False, "BMI_n": False, "Smoking_n": False, "Caloric_n": False}
        n_lines = 0    
        with open(To_open) as F:
                for line in F:
                        #if n_lines > 500: break
                        #n_lines += 1
                        l = line.rstrip().split("\t")
                        if "PSEUDOIDEXT" in line:
                                n = 0 
                                for Feature in l:
                                        if "PSEUDOIDEXT" in Feature: dic_entries["ID"] = n 
                                        if Feature == "AGE_1A1":  dic_entries["Age_n"] = n 
                                        elif Feature == "GESLACHT": dic_entries["Sex_n"] = n 
                                        elif Feature == "BMI": dic_entries["BMI_n"] = n 
                                        elif Feature == "SMK3": dic_entries["Smoking_n"] =  n
                                        elif Feature == "SUMOFKCAL": dic_entries["Caloric_n"] = n 
                                        n += 1
                        if l[dic_entries["ID"]] not in dic_covariates: dic_covariates[l[dic_entries["ID"]]] = {"Age_n":"NA","Sex_n":"NA", "BMI_n":"NA", "Smoking_n":"NA", "Caloric_n":"NA"}
                        for Feature in ["Age_n", "Sex_n", "BMI_n", "Smoking_n", "Caloric_n"]:
                                if dic_entries[Feature] != False:
                                        try: dic_covariates[l[dic_entries["ID"]]][Feature] = l[dic_entries[Feature]]
                                        except: continue
Header_cov = ["Age", "Sex", "BMI", "Smoking", "Caloric"]

l = len(dic_covariates)


Health_data = "/groups/umcg-lifelines/prm03/releases/pheno_lifelines/v1/tab_separated_labels/Questionnaire_Health_General.dat"

dic_sample = {}
n= 0
LINE = 0 

#HEALTH100E1
#HEALTH72C2
#HEALTH100E2
#HEALTH72C3

Questionnaire = [] #Questionnaire (Baseline, assesment 2, followup...)
Interest = ["HEALTH100E1","HEALTH72C2","HEALTH100E2","HEALTH72C3"]
with open(Health_data, "r") as File:
        for line in File:
                l = line.split("\t")
		#If first line, save the line as header
                if LINE == 0:  
                        header = l 
                        LINE = 1 
			interest_columns = []
			for n in range(4,len(header)):
				name = header[n]
				if name in Interest: interest_columns.append(n)
                        continue
                #if l[0] != "30387140": continue
                if "30" not in l[0]: continue
                if l[0] not in dic_sample: dic_sample[l[0]] = {}
		#Remove samples without all columns
                if len(l) != len(header):
                        continue
		#Go to each of the columns of interest
                for n in interest_columns: #range(4,len(header)):
                        column = l[n]
                        name = header[n]
                        #Create a dic of Sample: Time: Answers
                        if name not in dic_sample[l[0]]:
                                dic_sample[l[0]][name] = {}
                        column = column.strip()
                        if  not column: column = "NA"
   
                        if l[2] not in dic_sample[l[0]][name]:  
                                dic_sample[l[0]][name][l[2]] = column
                                #dic_sample[l[0]][name][l[2]].extend(to_add)        
                        if l[0] == "30387140" and l[2] not in Questionnaire: Questionnaire.append(l[2])
                        #print(column, name, l[0], l[2])
    

l = len(dic_sample)



HEAD = ["ID", "Disease", "Baseline assessment, questionnaire 1", "Second assessment, questionnaire 1", "Second assessment, questionnaire 2", "Follow-up questionnaire (1B)", "Follow-up questionnaire (1C)", "Age_n","Sex_n", "BMI_n", "Smoking_n"]

with open("TABLES/Data_v2.tsv", "w") as Final:
        Final.write("\t".join(HEAD) + "\n")
        for ID in dic_sample:
                for I in Interest:
                        Out = [ID]	
                        Out.append(I)
			#If the sample Lacks a Questionnaire
                        for n in Questionnaire:
                                try:  Out.append(dic_sample[ID][I][n])
                                except: Out.append("NaN")
                        to_add = []
                       
                        for keys in ["Age_n","Sex_n", "BMI_n", "Smoking_n"]:
                                try: ITEM = dic_covariates[ID][keys] 
                                except: print(ID) ; ITEM = "NA"
                                to_add.append(ITEM)                   
                        Out.extend(to_add)
                        
                        Entry = "\t".join(Out)
                        Final.write(Entry + "\n")


