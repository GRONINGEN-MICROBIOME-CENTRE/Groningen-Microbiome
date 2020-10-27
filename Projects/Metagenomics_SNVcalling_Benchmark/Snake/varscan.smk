
rule samtools:
	input:
		BAM=  Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		REF= lambda wildcards: Get_Ref(wildcards.ID)
	output:
		OUT = temp(Working_dir+"Raw_calls/{ID}_prevarscan.mpileup")
		
	shell:
		"ml SAMtools/1.9-foss-2015b; \n"
		"samtools mpileup -f {input.REF} {input.BAM} > {output.OUT} "


rule varscan:
	input: 
		rules.samtools.output.OUT
	output: 
		OUT = temp(Working_dir+"Raw_calls/{ID}_varscan.vcf")
	params:
		varscan = "/groups/umcg-lld/tmp04/umcg-haugustijn/SNV_calling//programs/VarScan2.jar"
	shell:
		"ml SAMtools/1.9-foss-2015b;\n"
		"java -jar {params.varscan} mpileup2snp {input} --output-vcf > {output.OUT} "		


#rule format_variants:
#	input: Working_dir+"Raw_calls/{ID}_varscan.vcf", #expand("{DIR}Raw_calls/{ID}_bcftools.vcf",DIR=Working_dir,ID = dic_ref.keys())
#	output: Working_dir+"Calls/varscan/{ID}.tsv" #expand('{DIR}Calls/{TOOL}/{ID}.tsv',DIR=Working_dir, TOOL = "bcftools",ID = dic_ref.keys())
#	shell:
#		"""
#		ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p varscan -o {output}
#		"""
