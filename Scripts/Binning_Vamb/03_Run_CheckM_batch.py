from pathlib import Path
from subprocess import call

Dic_all = {}
Overall = 0
counter = 0
for File in Path("vamb/bins/").glob("*.fna"):
	if counter >= 500:
		counter = 0
		Overall += 1 
	if Overall not in Dic_all: Dic_all[Overall] = []
	Dic_all[Overall].append(File)
	counter += 1

	
for key, Batch_files in Dic_all.items():
	key = str(key)
	if not Path("CheckM_batches/B_{key}".format(key=key)).exists():
		call("mkdir CheckM_batches/B_{key}".format(key=key), shell=True)
	for fna in Batch_files:
		N = fna.name
		NN = "CheckM_batches/B_{key}/".format(key=key) + N
		if Path(NN).exists(): continue
		Path(NN).symlink_to(fna.absolute())
		

for key, Batch_files in Dic_all.items():
	DIR =  "CheckM_batches/B_{key}".format(key=key)
	T = "5"
	O = "CheckM_batches/Completed/B_{key}".format(key=key)
	command = """#!/bin/sh
#SBATCH  --job-name=Check_M_batch_{M}.job
#SBATCH --output=logs/batch_{M}.out
#SBATCH --error=logs/batch_{M}.err
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH --cpus-per-task={T}
""".format(M="B"+str(key), T=T)
	command += "ml CheckM/1.1.2-foss-2019b-Python-3.7.4;\n" #"ml Anaconda3; source activate /data/umcg-tifn/Binning_vamb/Pipeline/Binning_env;\n"
	command += "checkm lineage_wf -f {output}/info -t {threads} -x fna {fna} {output}".format(output=O, threads=T, fna=DIR)
	with open("batch_script_{M}.sh".format(M=key), "w") as BATCH:
		BATCH.write(command)
	if key == "0": continue
	call("sbatch batch_script_{M}.sh".format(M=key), shell=True)
	call("rm batch_script_{M}.sh".format(M=key), shell=True)
