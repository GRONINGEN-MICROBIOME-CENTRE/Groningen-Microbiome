rule instrain:
	input:
		BAM=  Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		REF= lambda wildcards: Get_Ref(wildcards.ID)
	output:
		OUT = Working_dir+"Raw_calls/{ID}_instrain/output/{ID}_instrain_SNVs.tsv"
	params:
		instrain = "/groups/umcg-lld/tmp04/umcg-dwang/opt/Anaconda2/envs/py3/bin/inStrain",
		env= "/groups/umcg-lld/tmp04/umcg-dwang/opt/Anaconda2/envs/py3",
		Out_prefix =  Working_dir+"Raw_calls/{ID}_instrain"
	threads: 6
        shell:
		"ml Miniconda3/4.4.10 SAMtools/1.9-foss-2015b ; source activate {params.env} ;\n"
		"{params.instrain} profile {input.BAM}  {input.REF} -o {params.Out_prefix} -p {threads}"

#rule format_variants:
#        input: Working_dir+"Raw_calls/{ID}_instrain/output/{ID}_instrain_SNVs.tsv", #expand("{DIR}Raw_calls/{ID}_bcftools.vcf",DIR=Working_dir,ID = dic_ref.keys())
#        output: Working_dir+"Calls/instrain/{ID}.tsv" #expand('{DIR}Calls/{TOOL}/{ID}.tsv',DIR=Working_dir, TOOL = "bcftools",ID = dic_ref.keys())
#        shell:
#                """
#                ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p instrain -o {output}
#		"""


