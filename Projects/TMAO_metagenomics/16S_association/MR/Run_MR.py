dic_do = {}

with open("Associations.tsv") as F:
	for line in F:
		l = line.rstrip().split()
		dic_do[l[1]]= l[2]
		
from pathlib import Path
from subprocess import call
import sys
for Bug,metabolite in dic_do.items():
	Instrument_dir = "/groups/umcg-lld/tmp01/BOXY_BACKUP_umcg-lld/users/umcg-akurilshchikov/MiBioGen/Oct2020/betas/"
	Instrument = Instrument_dir+Bug+".summary.txt.gz"
	#Instrument = str(list(Instrument_dir.glob(Bug + ".summary.txt.gz"))[0])
	Outcome_dir = "Outcome/" #METAANALYSIS_CHOLINE_forSergio.txt.gz
	Outcome = Outcome_dir + "METAANALYSIS_"+metabolite.upper() + "_forSergio.txt.gz"
	#Outcome = str(list(Outcome_dir.glob("METAANALYSIS_"+metabolite.upper() + "_forSergio.txt.gz"))[0])
	
	script = "Run_{B}_{M}.sh".format(B=Bug, M=metabolite)
	Filter_intrument = "Instruments/{A}_raw.tsv".format(A=Bug)
	Work_dir = "Output_MR/{B}_{M}/".format(B=Bug, M=metabolite)
	if not Path(Work_dir).exists(): call("mkdir {W}".format(W=Work_dir), shell=True)
	with open(script, "w") as S:
		N = Bug+"_"+metabolite
		Text = """#!/bin/sh
#SBATCH  --job-name={N}.job
#SBATCH --output={N}.out
#SBATCH --error={N}.err
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=15G
#SBATCH --nodes=1\n""".format(N=N)
		if not Path(Filter_intrument).exists():
			Text += """Raw_instrumental={i}
zcat  $Raw_instrumental | awk 'NR==1; NR > 1 {{if ($10<5e-6) print $0}}' > {i2}
awk '{{print $0,"\\t",$2":"$3}}' {i2} >  /tmp/{bug}
mv  /tmp/{bug} {i2}\n""".format(i=Instrument,i2=Filter_intrument, bug=Bug)
		Text += "ml RPlus/3.6.1-foss-2018b-v20.10.1; Rscript  Perform_MR.R {I} {O} {D}".format(I=Filter_intrument, O=Outcome, D=Work_dir)
		S.write(Text)
	command = "sbatch "+script
	print(command)
	continue
	if len(sys.argv) > 1:
		if sys.argv[1] == "local": command = "bash "+script
	call(command, shell=True)
	
		



