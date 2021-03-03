#Idea:
#Use Alignment rules + calling rules + formatting rules from simulation pipeline
from pathlib import Path

#configfile:"config.yaml"

to_analyze = []
CHECK = "half_done.txt"
with open(CHECK) as F:
	for line in F:
		to_analyze.append(line.rstrip())


dic_fq = {}
for F in Path("FQ/NEW").glob("*.1.gz"):
	#SRS011134_kneaddata_merged.fastq.gz
	ID = F.name.split("_")[0]
	#if ID not in to_analyze: continue
	if ID not in dic_fq: dic_fq[ID] = []
	dic_fq[ID].append(str(F))
	dic_fq[ID].append(str(F).replace("1.gz","2.gz"))
	#break
Bacteria_to_analyze = []
with open("Bacterial_selection.tsv") as F:
	for line in F:
		Bacteria_to_analyze.append(line.split()[0])


dic_ref = {}
for Ref in Path("References").glob("*.fa"):
	ID = Ref.name.replace(".fa", "")
	if ID not in Bacteria_to_analyze: continue
	dic_ref[ID] = [str(Ref)]
def Get_Ref(ID):
	return(dic_ref[ID][0])
def Get_FQ(Sample, N):
	return(dic_fq[Sample][N])
def get_bacteria():
	return(dic_ref.values())
######COmmon rules for all samples (apply to the refernces·······

rule index_ref:
        input:
                lambda wildcards: Get_Ref(wildcards.ID)
        output:
                temp("References/indexed_{ID}")
        shell:
                "ml SAMtools/1.9-foss-2018b ;\n"
                "samtools faidx {input}"
rule Index_for_alignment:
	input:
                ref = lambda wildcards: Get_Ref(wildcards.ID)
	output:
		#flag = temp(touch("{ID}.indexed")),
		Index = "References/{ID}.1.bt2"
	params:
		save_dir = "References",
		outputpref = "{ID}"
	shell:
		"""
		module load Bowtie2/2.3.4.2-foss-2018b
		bowtie2-build {input.ref} {params.save_dir}/{params.outputpref}
		"""
#################################################################
#
#####Rules per sample: Alignment + VC###########################
#################################################################
#
#Alignment
rule Remove_human:
        input:
                fq1 = lambda wildcards: Get_FQ(wildcards.Sample, 0),
                fq2 = lambda wildcards: Get_FQ(wildcards.Sample, 1),
        output:
                fq1 = temp("{Sample}/Aligned/Cleaned_paired_1.fastq"),
                fq2 = temp("{Sample}/Aligned/Cleaned_paired_2.fastq"),
                log_file = temp("{Sample}/Logs/kneaddata")
        params:
                p = 7,
                outputdir = "{Sample}/Aligned/",
                db = "/groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/",
        threads: 7
        shell:
                "module load picard Python/2.7.11-foss-2015b Bowtie2/2.3.4.1-foss-2015b ;\n"
                "kneaddata -i {input.fq1} -i {input.fq2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix Cleaned --log {output.log_file}"

rule Merge_ref:
	input: get_bacteria()
	output: 
		ref = "Uniref.fa",
		Index = "Uniref.1.bt2"
	shell:
		"cat {input} > {output.ref} ;\n"
		"module load Bowtie2/2.3.4.2-foss-2018b ;\n"
		"bowtie2-build {output.ref} Uniref"
rule Alignment:
        input:
                ref = rules.Merge_ref.output.ref,
                F1 = lambda wildcards: Get_FQ(wildcards.Sample, 0),
                F2 = lambda wildcards: Get_FQ(wildcards.Sample, 1),
        output:
                sam = temp("{Sample}/Aligned/temp.sam"),
                bam = temp("{Sample}/Aligned/temp.bam"),
                sorted_bam = temp("{Sample}/Aligned/Sorted.bam")
        threads: 12
        shell:
                """
                module load Bowtie2/2.3.4.2-foss-2018b
                bowtie2 --rg-id {wildcards.Sample} --rg SM:{wildcards.Sample} --rg LIB:{wildcards.Sample} --rg PL:Illumina -p {threads} -x Uniref -1 {input.F1} -2 {input.F2} -S {output.sam}

                module load picard SAMtools/1.9-foss-2018b
                samtools view -bS {output.sam} > {output.bam}
                samtools sort {output.bam} > {output.sorted_bam}
                """

rule add_group:
        input:
                BAM = rules.Alignment.output.sorted_bam
        output:
                BAM =  temp("{Sample}/Aligned/temp_clean2.bam"),
        shell:
                "ml picard SAMtools/1.9-foss-2018b ;\n"
                "java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={input.BAM} O={output.BAM} RGID={wildcards.Sample} RGLB={wildcards.Sample} RGPL=ILLUMINA RGPU=unit1 RGSM={wildcards.Sample} ;\n"
                

rule MarkDuplicates:
        input:
                bam = rules.add_group.output.BAM #"{Sample}/Aligned/Sorted.bam",
        output:
                derep  = temp("{Sample}/Aligned/clean_temp.bam"),
                metrics = temp("{Sample}/Logs/removedup"),
                BAM = "{Sample}/Aligned/clean.bam",
		BAI =  "{Sample}/Aligned/clean.bam.bai"
        params:
                tempdir = "{Sample}/tmp/",
                Picath = "/apps/software/picard/2.20.5-Java-11-LTS"
        shell:
                """
                module load picard
                java -jar {params.Picath}/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output.derep} REMOVE_DUPLICATES=True METRICS_FILE={output.metrics}
                java -jar {params.Picath}/picard.jar CleanSam INPUT={output.derep} OUTPUT={output.BAM} TMP_DIR={params.tempdir}
		ml SAMtools/1.9-foss-2018b
		samtools index  -b {output.BAM}
		
                """




#############################
### Variatn Calling ########################
def Variant_input(ID, tool):
	if tool in ["bcftools", "freebayes","mutect2", "haplotypecaller","varscan"]:
		INPUT = "{ID}/call_{tool}.vcf".format(ID=ID, tool = tool)
	elif tool == "instrain":
		INPUT =  "{Sample}/{Sample}_instrain/output/{Sample}_instrain_SNVs.tsv".format(Sample=ID)
	return(INPUT)

include: "Snake/gatk.smk"
include: "Snake/instrain.smk"

rule Format:
	input: lambda wildcards: Variant_input(wildcards.Sample, wildcards.tool)
	output: "{Sample}/Calls/{tool}.tsv"
	shell:
		"python3 Reformat.py  -f {input} -p {wildcards.tool} -o {output}"
	



calls = [expand("{Sample}/{Sample}_instrain/output/{Sample}_instrain_SNVs.tsv", Sample= dic_fq.keys()), expand("{Sample}/call_mutect2.vcf", Sample= dic_fq.keys()), expand("{Sample}/call_haplotypecaller.vcf", Sample=  dic_fq.keys())]



rule all:
	input:
		#index = expand("{ID}.indexed", ID=dic_ref.keys())
		Alignment = expand("{Sample}/Aligned/clean.bam", ID = dic_ref.keys(), Sample = dic_fq.keys()),
		#Raw_call = calls,
		Calls = expand("{Sample}/Calls/{tool}.tsv",  Sample = dic_fq.keys(), tool=["mutect2","haplotypecaller", "instrain"])
		#calls = expand('{Sample}/Calls/{TOOL}/{ID}.tsv', TOOL = Tools, Sample = dic_ref.keys(), ID =Bacteria_taxa)



