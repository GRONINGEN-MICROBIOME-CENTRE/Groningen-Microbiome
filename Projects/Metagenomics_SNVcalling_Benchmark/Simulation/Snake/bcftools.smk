
rule bcftools:
	input:
		BAM=  Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		REF= lambda wildcards: Get_Ref(wildcards.ID)
	output:
		OUT = Working_dir+"Raw_calls/{ID}_bcftools.vcf"
		
	shell:
		"ml SAMtools/1.9-foss-2015b BCFtools/1.7-foss-2015b ;"
		"bcftools mpileup -f {input.REF} {input.BAM} | bcftools call --ploidy 1 -mv  -o {output.OUT} "
#rule format_variants:
#	input: Working_dir+"Raw_calls/{ID}_bcftools.vcf", #expand("{DIR}Raw_calls/{ID}_bcftools.vcf",DIR=Working_dir,ID = dic_ref.keys())
#	output: Working_dir+"Calls/bcftools/{ID}.tsv" #expand('{DIR}Calls/{TOOL}/{ID}.tsv',DIR=Working_dir, TOOL = "bcftools",ID = dic_ref.keys())
#	shell:
#		"""
#		ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p bcftools -o {output}
#		"""
