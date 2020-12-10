
rule freebayes:
	input:
		BAM=  Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		REF= lambda wildcards: Get_Ref(wildcards.ID)
	output:
		OUT_temp = temp(Working_dir+"Raw_calls/{ID}_freebayes_temp.vcf"),
		OUT = Working_dir+"Raw_calls/{ID}_freebayes.vcf"
		
	shell:
		"""
		ml freebayes/1.1.0-foss-2015b-961e5f3
		freebayes -f {input.REF} --ploidy 1 --min-alternate-count 2 --min-alternate-fraction 0 --pooled-continuous {input.BAM} > {output.OUT_temp}
		ml Miniconda3/4.4.10 ; source activate /groups/umcg-lld/tmp04/umcg-sandreusanchez/Conda_env/METASNV
		vt decompose_blocksub -a -o {output.OUT} {output.OUT_temp}
		"""
#rule format_variants_freebayes:
#	input: Working_dir+"Raw_calls/{ID}_freebayes.vcf",
#	output: Working_dir+"Calls/freebayes/{ID}.tsv"
#	shell:
#		"""
#		 ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p freebayes -o {output}
#		"""
