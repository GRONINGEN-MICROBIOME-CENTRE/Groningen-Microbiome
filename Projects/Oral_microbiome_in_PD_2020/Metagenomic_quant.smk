from pathlib import Path
from subprocess import call
basedir = workflow.basedir
configfile: "config.yaml"

###Configuration should have a path to a file where each line represents a fastq file to be analyzed (both forward and reverse per line)
F = config["Paths_file"]
Taxonomy = config["Taxonomy"]




dic_files = {}


with open(F) as Interest_paths:
	for path in Interest_paths:
		P = path.split(",")
		ID = P[0]
		if not Path("Results/"+ID).exists():
			call("mkdir Results/"+ID,shell=True)
		if not Path("tmp/"+ID).exists():
			call("mkdir tmp/"+ID,shell=True)
		Forward = P[1]
		Reverse = P[2].rstrip()
		dic_files[ID]= [Forward, Reverse]
		### Set up a naming format, so that I can get an ID name that shall be related with both forward and reverse reads

### Set up output: Matrix of counts (Indv one, we might add one more step for merging)
localrules : all, Merge_reads, Merge_Tables, mv_pathways, Merge_taxonomy

rule all:
	input:
		expand('Results/{ID}/{ID}_feature_counts.tsv',zip,ID=dic_files.keys()),
		expand('Results/{ID}/Done.flag',zip,ID=dic_files.keys()),
		"Results/Merged_batch_pathabundance.tsv",
		"Results/Merged_batch.tsv"

def find_reads(ID, dire):
	D = dic_files[ID]
	if dire == "f":
		return(D[0])
	else:
		return(D[1])


	
	
### Steps are based on Arnau Vich's pipeline in https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome/blob/master/Scripts/Metagenomics_pipeline_v1.md

rule RemoveHuman:
	input:
		forward = lambda wildcards: find_reads(wildcards.ID, dire="f"),
		reverse = lambda wildcards: find_reads(wildcards.ID, dire="r")
	output:
		clean_f =  "tmp/{ID}/filtering_Data/{ID}_paired_1.fastq",
		clean_r = "tmp/{ID}/filtering_Data/{ID}_paired_2.fastq"
	params:
		db = config["Homo_sapiens_DB"] #"/groups/umcg-gastrocol/tmp04/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/",
		nthreads = 6,
		nprocesses = 7,
		trimmomatic =  config["trimmomatic"] #"/home/umcg-sandreusanchez/.conda/envs/myenv/share/trimmomatic"
	log:
		"logs/kneaddate_{ID}"
	shell:
		"""
		module load Trimmomatic/0.32-Java-1.7.0_80
		module load picard 
		module load Python/3.4.1-foss-2015b
		module load Bowtie2
		module load kneaddata
		kneaddata --input {input.forward} --input {input.reverse} -db {params.db} --output tmp/{wildcards.ID}/filtering_Data/ --log {log} -t {params.nthreads} -p {params.nprocesses} --trimmomatic {params.trimmomatic} --output-prefix {wildcards.ID}
		"""
rule Merge_reads:
	input:
		f = rules.RemoveHuman.output.clean_f,
		r = rules.RemoveHuman.output.clean_r
	output:
		merged = "tmp/{ID}/{ID}_merged.fq"
	shell:
		"cat  {input.f} {input.r} > {output.merged}"

rule Taxon_Classification:
	input:
		Reads = rules.Merge_reads.output.merged
	output: 	
		'Results/{ID}/{ID}_feature_counts.tsv'
	params:
		metadata = config["Metaphlan_metadata] "#"/groups/umcg-gastrocol/tmp04/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200.pkl", #?
		process = 6,
		metaphlan2 = config["metaphlan"] #"/groups/umcg-gastrocol/tmp04/metagenomic_tools/metaphlan_2/metaphlan2.py" 
	shell:
		"""
		module load Python/3.4.1-foss-2015b
		mkdir tmp/{wildcards.ID}/metaphlan
		module load Bowtie2
		{params.metaphlan2} {input.Reads}  --input_type multifastq --mpa_pkl {params.metadata} --nproc {params.process} -o {output} --tmp_dir tmp/{wildcards.ID}/metaphlan
		"""


rule Merge_Tables:
	input: expand('Results/{ID}/{ID}_feature_counts.tsv',zip,ID=dic_files.keys())
		
	output:	"Results/Merged_batch.tsv"
	
	params:
		Merge = "/groups/umcg-gastrocol/tmp04/metagenomic_tools/metaphlan_2/utils/merge_metaphlan_tables.py"
	shell:
		"{params.Merge} {input} > {output}"

if Taxonomy == True:
	rule functional_profiling:
		input:
			Reads = rules.Merge_reads.output.merged,
			taxonomy = rules.Taxon_Classification.output
		output:
			directory = directory('Results/{ID}/{ID}_pathway_counts'),
			flag = touch('Results/{ID}/Done.flag')
		params:
			process = 6,
			metaphlan2 =  config["metaphlan"] #"/groups/umcg-gastrocol/tmp04/metagenomic_tools/metaphlan_2/metaphlan2.py",
			tlog = "tmp/{ID}/{ID}.full.humann2.log"
		shell:
			"""
			module load picard
			module load Python/3.4.1-foss-2015b
			module load Bowtie2
			ml humann2
			humann2 --input {input.Reads} --output {output.directory} --taxonomic-profile {input.taxonomy} --threads {params.process} --o-log {params.tlog} --remove-temp-output
			"""

	rule mv_pathways:
		input:
			rules.functional_profiling.output.directory
		output:
			'Results/Pathway_tables/{ID}.pathabundance.tsv'
		shell:
			"""
			[ ! -d  "Results/Pathway_tables/" ] && mkdir Results/Pathway_tables/
			cp {input}/{wildcards.ID}_merged_pathabundance.tsv {output}
			"""
	rule Merge_taxonomy:
		input:  expand('Results/Pathway_tables/{ID}.pathabundance.tsv',zip,ID=dic_files.keys())
		output:	 "Results/Merged_batch_pathabundance.tsv"
		shell:
			"""
			ml humann2	
			humann2_join_tables --input Results/Pathway_tables/ --output {output}
			"""
