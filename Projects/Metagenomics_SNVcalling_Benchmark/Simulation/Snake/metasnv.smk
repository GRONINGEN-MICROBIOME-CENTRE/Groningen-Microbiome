
#export PATH=\$PATH:/home/umcg-lchen/bin/
#export PATH=\$PATH:/groups/umcg-lld/tmp04/umcg-lchen/metaSNV/metaSNV_install/




rule metasnv:
	input:
		BAM=  Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		REF= lambda wildcards: Get_Ref(wildcards.ID)
	output:
		BAM_path = temp(Working_dir+"Raw_calls/{ID}_metasnv_path.txt"),
		#OUT = Working_dir+"Raw_calls/MetaSNV/{ID}/{ID}.all_cov.tab",
		#Call = directory(Working_dir+"Raw_calls/MetaSNV/{ID}/filtered/filtered/pop"),
		Call = Working_dir + "Raw_calls/MetaSNV/{ID}/snpCaller/called_SNPs"
	params:
		metaSNV = "/groups/umcg-lld/tmp04/umcg-sandreusanchez/Conda_env/METASNV/bin/metaSNV.py", #"/groups/umcg-lld/tmp04/umcg-lchen/metaSNV/metaSNV_install/metaSNV.py",
		metaSNV_post = "/groups/umcg-lld/tmp04/umcg-lchen/metaSNV/metaSNV_install/metaSNV_post.py",
		Directory = Working_dir+"Raw_calls/MetaSNV/{ID}"
	shell:
		"""
		ml Miniconda3/4.4.10
		source activate /groups/umcg-lld/tmp04/umcg-sandreusanchez/Conda_env/METASNV
		
		if [ -d "{params.Directory}" ]; then rm -Rf {params.Directory}; fi
	

		ml ncurses/5.9-foss-2015b Python/2.7.11-foss-2015b pandas/0.18.1-foss-2015b-Python-2.7.11 numpy/1.11.0-foss-2015b-Python-2.7.11 SAMtools/0.1.19-foss-2015b HTSlib/1.9-foss-2015b
		
		echo {input.BAM} > {output.BAM_path}
	
		metaSNV.py  {params.Directory} {output.BAM_path} {input.REF}
		
		"""
			
#rule format_variants_metasnv:
#	input: rules.metasnv.output.Call
#	output: Working_dir+"Calls/metasnv/{ID}.tsv"
#	shell:
#		"ml Python/3.4.1-foss-2015b ; python3 scripts/Reformat.py -f {input} -p freebayes -o {output}"
