#/groups/umcg-lld/tmp04/umcg-dwang/SNV/01.alignment/
def Get_Fasq(read=1):
	if config["fastq"] == "NA":
		if read == 1:  Fq = Working_dir + "Simulation/FQ/simulated_R1.fastq"
		elif read == 2: Fq = Working_dir + "Simulation/FQ/simulated_R2.fastq"
	else:
		Fqs = config["fastq"].split(",")
		Fq = Fqs[read-1]
	return(Fq)

for ID in dic_ref:
	DIR = Path(Working_dir+"Alignment/"+ID)
	if not DIR.exists(): DIR.mkdir(parents=True, exist_ok=True)
rule Index:
	input: 
		ref = lambda wildcards: Get_Ref(wildcards.ID),
	output: 
		flag = temp(touch(Working_dir+"{ID}.indexed")),
		Index = Working_dir+"Alignment/{ID}/{ID}.1.bt2"
	params:
		save_dir = Working_dir+"Alignment/{ID}",
		outputpref = "{ID}"
	shell:
		"""
		module load Bowtie2/2.3.4.1-foss-2015b
		if [ ! -d "{params.save_dir}" ] ; then mkdir {params.save_dir} ; fi
		sleep 5
		bowtie2-build {input.ref} {params.save_dir}/{params.outputpref}
		"""
rule Remove_human:
	input: 
		fq1 = Get_Fasq(1),
		fq2 = Get_Fasq(2),
		ready_flag = Working_dir + "Simulations_completed.flag"
	output:
		fq1 = Working_dir+"Alignment/Cleaned_paired_1.fastq",
		fq2 = Working_dir+"Alignment/Cleaned_paired_2.fastq",
		log_file = Working_dir+"Logs/kneaddata"
	params:
		p = 7,
		outputdir = Working_dir+"Alignment/",
		db = "/groups/umcg-lld/tmp04/umcg-dwang/database/keaddata-0.5.4",
	threads: 7
	shell:
		"module load picard Python/2.7.11-foss-2015b Bowtie2/2.3.4.1-foss-2015b ;\n" 
		"kneaddata -i {input.fq1} -i {input.fq2} -db {params.db} -t {threads} -p {params.p} -o {params.outputdir} --output-prefix Cleaned --log {output.log_file}" 

rule Alignment:
	input: 
		ref = rules.Index.output.flag,
		F1 = rules.Remove_human.output.fq1,
		F2 = rules.Remove_human.output.fq2,
	output:
		sam = temp(Working_dir+"Alignment/Aligned/temp_{ID}.sam"),
		bam = temp(Working_dir+"Alignment/Aligned/temp_{ID}_unclean.bam"),
		sorted_bam = temp(Working_dir+"Alignment/Aligned/temp_{ID}.bam")
	params:
		index = Working_dir+"Alignment/{ID}"
	threads: 12
	shell:
		"""
  		module load Bowtie2/2.3.4.1-foss-2015b
		bowtie2 --rg SM:mock_community --rg LIB:{wildcards.ID} --rg PL:Illumina -p {threads} -x {Working_dir}Alignment/{wildcards.ID}/{wildcards.ID} -1 {input.F1} -2 {input.F2} -S {output.sam}

                module load picard SAMtools/1.9-foss-2015b
                samtools view -bS {output.sam} > {output.bam}
                samtools sort {output.bam} > {output.sorted_bam}
		"""
rule MarkDuplicates:
	input:
		bam = rules.Alignment.output.sorted_bam,
	output:
		derep  = temp(Working_dir+"Alignment/Aligned/{ID}_clean_temp.bam"),
		metrics = Working_dir+"Logs/removedup_{ID}",
		derep2 = Working_dir+"Alignment/Aligned/{ID}_clean.bam"
	params:
		tempdir = Working_dir+"tmp/",
		Picath = "/apps/software/picard/2.18.26-Java-1.8.0_74"
	shell:
		"""
		module load picard
                java -jar {params.Picath}/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output.derep} REMOVE_DUPLICATES=True METRICS_FILE={output.metrics}
                java -jar {params.Picath}/picard.jar CleanSam INPUT={output.derep} OUTPUT={output.derep2} TMP_DIR={params.tempdir}
		"""	
