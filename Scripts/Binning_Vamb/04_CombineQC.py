from pathlib import Path
import pandas as pd
import re

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Bin Id                 Marker lineage             # genomes   # markers   # marker sets    0     1     2    3   4   5+   Completeness   Contamination   Strain heterogeneity  
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  S89C58        c__Deltaproteobacteria (UID3218)        61        
Complete = 50
Contamination = 5
QC_s = 50
#Params taken from https://www.nature.com/articles/s41587-020-0603-3#Sec11

Header = ["Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "0",  "1" , "2", "3", "4", "5+", "Completeness", "Contamination", "Strain heterogeneity", "QC_pass"]
Table = None
for F in Path("CheckM_batches/Completed/").glob("*/info"):
	with open(F) as lines:
		for line in lines:
			if line[0:2] == "--": continue
			if "Bin Id" in line:
				continue
			else:
				Fields  = re.split(r'\s{2,}', line.strip())
				Com = float(Fields[-3]) ; Cont = float(Fields[-2]) ; Q= Com-5*Cont
				if Com > Complete and Cont < Contamination and  Q > QC_s: Pass =  True
				else: Pass=False
				Fields.append(Pass)
				Fields = pd.DataFrame(Fields).T
				Fields.columns = Header
				if Table is None:
					Table = Fields
				else:
					Table = Table.append(Fields)
Table.to_csv("vamb/QC_check.csv")
