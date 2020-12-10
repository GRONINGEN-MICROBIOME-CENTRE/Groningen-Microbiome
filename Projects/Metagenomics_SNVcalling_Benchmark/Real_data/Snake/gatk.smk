rule prepare_Dic:
	input:
		Ref = "Uniref.fa"	
	output:
		fai =  "Uniref.fa.fai",
		dicti = "Uniref.dict",
	shell:
		"ml picard/2.20.5-Java-11-LTS SAMtools/1.9-GCCcore-7.3.0 ;\n"
		"samtools faidx {input.Ref} -o {output.fai} ;\n"
		"java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R={input.Ref} O={output.dicti}"

rule mutect2:
	input: 
		BAM=rules.MarkDuplicates.output.BAM,
		BAI=rules.MarkDuplicates.output.BAI,
		REF="Uniref.fa",
		FAI="Uniref.fa.fai"
	output:
		mutect = "{Sample}/call_mutect2.vcf"
	shell:
		"""
		ml GATK/4.1.4.1-Java-8-LTS
		gatk Mutect2 --af-of-alleles-not-in-resource 0.33 -R {input.REF} -I {input.BAM} -O {output.mutect}
		"""
rule haplotypecaller:
	input:
		BAM=rules.MarkDuplicates.output.BAM,
		BAI=rules.MarkDuplicates.output.BAI,
		REF="Uniref.fa",
		FAI="Uniref.fa.fai"
	output:
		hapcall = "{Sample}/call_haplotypecaller.vcf"
		
	shell:
		"""
		ml GATK/4.1.4.1-Java-8-LTS
		ml picard/2.20.5-Java-11-LTS SAMtools/1.9-GCCcore-7.3.0
		gatk HaplotypeCaller -ploidy=1 -R {input.REF} -I {input.BAM} -O {output.hapcall}
		"""
