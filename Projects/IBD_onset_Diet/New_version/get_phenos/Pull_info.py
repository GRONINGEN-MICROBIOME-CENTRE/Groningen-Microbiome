import sys

#1.  Get input 
try:
	Item = sys.argv[1]
	if len(sys.argv) != 3: Output = "Output_script.tsv"
	else: Output = sys.argv[2]
except: print("To use the script please provide the input and output files")


#2. Generate a dic with the files to open and the fields to look for

Items = {}
with open(Item) as F:
	#grep file
	for line in F:
		#Take only lines with :
		if not ":" in line: continue
		line = line.rstrip()
		l = line.split(":")
		if len(l) < 2: continue
		#Get the field name
		Name = l[1].split(",")[0]
		if len(Name) < 2: continue
		#Get the file name
		File = l[0]
		if File not in Items:
			Items[File]  = []
		#Save fileds in file name
		Items[File].append(Name)

# Go participant by participant in each of the files in Items, and take the fields named as each of the items in the dic
dic_participant = {}
Header = ["ID"]
Times = {}

for File in Items:
	print("Iteration in File: " + File)
	#Getting info to call the column regarding the time of the followup, and if its a follow-up quest.
	Time = File.split("/")[-1][0]
	Quest_n = File.split("/")[-1][1]
	#if Quest_n != "b": continue
	if Quest_n != "a": Time = Time + "." + Quest_n
	#In each timetable add the Questions of interest
	Times[Time] = []
	with open(File) as F:
		#Counter to know when is the header
		n = 0
		list_items = []
		item_name = Items[File]
		print("Items in search file: "+ " ".join(item_name))
		for line in F:
			l = line.rstrip().split(",")
			if n == 0 :
				# IF header, save the line so we have a reference of how long should each row be
				H = l
				for item in item_name:
					#Add item name together with time in the header line
					Header.append(item.strip('"') + "_"+Time)
					#Search for the index of a field
					list_items.append(l.index(item))
					Times[Time].append(item)
				n +=1
				print(list_items)
			else:
				participant = l[0]
				if participant != '"e3d29040-40d8-4f58-b66a-32769cdba2bd"': continue
				if participant not in dic_participant: dic_participant[participant] = {}
				if Time not in dic_participant[participant]: dic_participant[participant][Time] = []
				#If length of the row is not good, this samples is not reliable, add NAs
				print(len(l), len(H))
				if len(l) != len(H):
					for i in list_items:
						dic_participant[participant][Time].append("NA")
					continue
				#if participant != '"380c7a43-33bf-4357-84b1-6d8d3e007233"': continue
				#If length was good
				for i in list_items:
					# i is the index of the items we look for 
					Add = l[i].strip('"')
					#print(Add)
					if "$" in Add: Add = "NA"
					dic_participant[participant][Time].append(Add)
#print(dic_participant['"380c7a43-33bf-4357-84b1-6d8d3e007233"'])


exit()
#Times_h = []
#['ID', 'colitisulcerosa_followup_adu_q_1_1.b', 'crohns_followup_adu_q_1_1.b', 'crohns_presence_adu_q_1_1', 'colitisulcerosa_presence_adu_q_1_1', 'colitisulcerosa_followup_adu_q_1_2', 'crohns_followup_adu_q_1_2', 'colitisulcerosa_followup_adu_q_1_1.c', 'crohns_followup_adu_q_1_1.c', 'colitisulcerosa_followup_adu_q_1_3', 'crohns_followup_adu_q_1_3']
#for i in Header:
#	if i == "ID": continue
#	Time = 
#	Times_h = 

with open(Output, "w") as O:
	H = "\t".join(Header) + "\n"
	O.write(H)
	for participant  in dic_participant:
		#For each participant, for each time
		Q = []
		for Time in Times:
			Expected = len(Times[Time]) #Number of Questions of interst in that particular time
			Ques = dic_participant[participant]
			if Time in Ques: Ques = Ques[Time] #If participnat had question in that time
			else: #Othersie, fill with NAs
				Ques = []
				[Ques.append("NA") for i in range(Expected)]
			Q.extend(Ques)
		if participant == '"380c7a43-33bf-4357-84b1-6d8d3e007233"': print(Header,Q)
		W = participant + "\t" + "\t".join(Q) + "\n"
		O.write(W)








