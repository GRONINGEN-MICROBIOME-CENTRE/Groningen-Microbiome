from pathlib import Path
from subprocess import call

#snakemake -s Run_shortbred_process.smk -np all

basedir = workflow.basedir
configfile: "config.yaml"
localrules: All_stats_together, Index
Output_list = []

dic_ref = {}

#Samples, format: Ref
Ref = config["Ref_list"]
#Working dir were everything will be stored
Working_dir = config["Directory"]

if Working_dir[-1] != "/": Working_dir = Working_dir + "/"

if not Path(Working_dir).exists():
	Path(Working_dir).mkdir(parents=True, exist_ok=True)

if not Path(Ref).exists():
	print("Error in config.yaml: 'Ref_list': {SAMPLE} was not found".format(SAMPLE=Ref))
	exit()
#Do simulations
Do_simulation = config["simulation"]
if Do_simulation not in [True, False]:
	print("Error in config.yaml: 'simulation' value should be set to either 'True' or 'False' (with no quote marks)")
	exit()

#Create ID dic:
with open(Ref) as F:
	for line in F:
		Ref = line.rstrip()
		ID = Path(line.rstrip()).stem
		dic_ref[ID] = [Ref]
Working_dir+"Alignment/Aligned/{ID}_clean.bam"
if Do_simulation == True:
	Output_list.append("{PATH}Simulations_completed.flag".format(PATH=Working_dir))
else:
	with open("{PATH}Simulations_completed.flag".format(PATH=Working_dir), "w") as T: pass 
	if "," in config["fastq"]: 
		print("Error in config.yaml: 'fastq' value should include two paths for paired-end fastq separated by a ','.")
		exit()
design = config["design"]
#Check tools to be used
Available_tool = ["all" ,"bcftools","mutect2", "haplotypecaller", "freebayes", "instrain", "metasnv", "varscan"]
tools = config["tools"]
if tools == "all": Tools = Available_tool[1:]
else: Tools = tools.split(",")
for tool in tools.split(","):
	if tool not in Available_tool:
		print("Error: Tool {NAME} provided in 'tools' is not in the list of available tools {LIST}".format(NAME=tool, LIST=",".join(Available_tool)))
		exit()
			
##TARGETS##
rule all:
	input:	
		Output_list,
		Alignment = expand(Working_dir+"Alignment/Aligned/{ID}_clean.bam", ID = dic_ref.keys()),
		calls = expand('{DIR}Calls/{TOOL}/{ID}.tsv',DIR=Working_dir, TOOL = Tools, ID = dic_ref.keys()),
		Stats_calls = expand("{DIR}Statistics_calling/{TOOL}.{ID}.tsv",DIR=Working_dir, TOOL = Tools, ID = dic_ref.keys()),
		ROC_curves = expand("{DIR}Statistics_calling/{TOOL}.{ID}.tsvROC",DIR=Working_dir, TOOL = ["haplotypecaller","freebayes","bcftools"], ID = dic_ref.keys()),
		Done_flag = Working_dir+"Pipeline_ended.flag"
if Do_simulation == True:
	include: "Snake/Simulation.smk"
include: "Snake/Alignment.smk"

rule index_ref:
	input:
		lambda wildcards: Get_Ref(wildcards.ID)
	output:
		temp(touch(Working_dir+"tmp/indexed_{ID}"))
	shell:
		"ml SAMtools/1.9-foss-2015b ;\n"
		"samtools faidx {input}"

for tool in Tools:
	if tool == "mutect2" or tool == "haplotypecaller": tool = "gatk"
	Script = "Snake/{tool}.smk".format(tool=tool)
	include: Script

def Variant_input(ID, tool):
	if tool in ["bcftools", "freebayes","mutect2", "haplotypecaller","varscan"]:
		INPUT = Working_dir+"Raw_calls/{ID}_{tool}.vcf".format(ID=ID, tool = tool)
	elif tool == "instrain":
		INPUT = Working_dir+"Raw_calls/{ID}_instrain/output/{ID}_instrain_SNVs.tsv".format(ID=ID)
	elif tool == "metasnv":
		INPUT =  Working_dir + "Raw_calls/MetaSNV/{ID}/snpCaller/called_SNPs".format(ID= ID)
	return(INPUT)
rule format_variants:
        input: lambda wildcards: Variant_input(wildcards.ID, wildcards.tool) #Working_dir+"Raw_calls/{ID}_bcftools.vcf",
        output: Working_dir+"Calls/{tool}/{ID}.tsv"
        shell:
                """
                ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p {wildcards.tool} -o {output}
                """
include: "Snake/Check_Calls.smk"

rule All_stats_together:
	input: 
		Stats = expand("{DIR}Statistics_calling/{TOOL}.{ID}.tsv",DIR=Working_dir, TOOL = Tools, ID = dic_ref.keys()),
	output:
		Stats_together = temp(Working_dir+"tmp/Stats_together.tsv"),
		Benchmark = Working_dir+"Benchmark.tsv",
		Done_flag = touch(Working_dir+"Pipeline_ended.flag")
	run:
		shell("cat {input.Stats} > {output.Stats_together}")
		with open(output.Benchmark, "w") as Out: pass
		dic_ref2 = {}
		for ID  in dic_ref:
			total_size = []
			Ref = dic_ref[ID]
			for Item in Ref:
				size = 0
				with open(Item) as R:
					for line in R:
						if line[0] == ">": continue
						size += len(line.rstrip())
				total_size.append(str(size))
			dic_ref2[ID] = ",".join(total_size)
		with open(output.Stats_together) as I:
			m = 0
			for line in I:
				if m == 0: line = line.rstrip() +"\t" +"\t".join(["Length_genome", "Benchmark_Design", "\n"])
				else:
					if "TP" in line: continue
					ID = line.rstrip().split()[-2]
					Length_dic = dic_ref2[ID]
					Number_ref = design
					line = line.rstrip() +"\t" +"\t".join([Length_dic, Number_ref, "\n"])
				m += 1
				with open(output.Benchmark, "a") as Out:
					Out.write(line)
					
