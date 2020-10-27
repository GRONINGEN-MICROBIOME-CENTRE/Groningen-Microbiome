rule prepare_Dic:
	input:
		REF = lambda wildcards: Get_Ref(wildcards.ID),
		Fai = Working_dir+"tmp/indexed_{ID}"
	output:
		ref2 = Working_dir+"Raw_calls/REF/{ID}.fa",
		dicti = Working_dir+"Raw_calls/REF/{ID}.dict",
	shell:
		"ml picard/2.18.26-Java-1.8.0_74 ;\n"
		"PATH_i=$(readlink -f  {input.REF}) ;\n"
		"ln -s $PATH_i {output.ref2} ;\n"
		"ln -s  $PATH_i.fai  {output.ref2}.fai ;\n"
		"java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R={output.ref2} O={output.dicti}"

rule add_group:
	input:
		BAM = Working_dir+"Alignment/Aligned/{ID}_clean.bam"
	output:
		BAM =  temp(Working_dir+"Alignment/Aligned/{ID}_clean2.bam"),
		BAI = temp(Working_dir+"Alignment/Aligned/{ID}_clean2.bam.bai"),
	shell:
		"ml picard/2.18.26-Java-1.8.0_74 SAMtools/1.9-foss-2015b ;\n"
		"java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I={input.BAM} O={output.BAM} RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20 ;\n"
		"samtools index {output.BAM} {output.BAI}"

rule gatk:
	input:
		BAM=  rules.add_group.output.BAM, #Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		BAI = rules.add_group.output.BAI,
		REF= rules.prepare_Dic.output.ref2,
	output:
		mutect = Working_dir+"Raw_calls/{ID}_mutect2.vcf",
		hapcall = Working_dir+"Raw_calls/{ID}_haplotypecaller.vcf"
		
	shell:
		"""
		ml Miniconda3/4.4.10
		source activate /groups/umcg-lld/tmp04/umcg-sandreusanchez/Conda_env/METASNV
		ml picard/2.18.26-Java-1.8.0_74 SAMtools/1.9-foss-2015b
		gatk HaplotypeCaller -ploidy=1 -R {input.REF} -I {input.BAM} -O {output.hapcall}
		gatk Mutect2 --af-of-alleles-not-in-resource 0.33 -R {input.REF} -I {input.BAM} -O {output.mutect}
		"""
			
#rule format_variants_mutect2:
#	input: Working_dir+"Raw_calls/{ID}_mutect2.vcf",
#	output: Working_dir+"Calls/mutect2/{ID}.tsv"
#	shell:
#		" ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p mutect2 -o {output}"
#rule format_variants_haplocaller:
#	input: Working_dir+"Raw_calls/{ID}_haplotypecaller.vcf"
#	output: Working_dir+"Calls/haplotypecaller/{ID}.tsv"
#	shell:
#		" ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p haplotypecaller -o {output}"
