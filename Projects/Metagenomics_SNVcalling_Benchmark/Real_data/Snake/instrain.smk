rule instrain:
	input:
		BAM=rules.MarkDuplicates.output.BAM,
		REF="Uniref.fa",
	output:
		OUT ="{Sample}/{Sample}_instrain/output/{Sample}_instrain_SNVs.tsv"
	params:
		instrain = "/groups/umcg-tifn/tmp01/users/umcg-dwang/opt/Anaconda3/envs/instrain/bin/inStrain",
		env=  "/groups/umcg-tifn/tmp01/users/umcg-dwang/opt/Anaconda3/envs/instrain",
		Out_prefix =  "{Sample}/{Sample}_instrain"

	threads: 6
        shell:
		"ml Anaconda3/5.3.0 SAMtools/1.9-GCCcore-7.3.0 ; source activate {params.env} ;\n"
		"""
		{params.instrain} profile {input.BAM}  {input.REF} -o {params.Out_prefix} -p {threads}
	
		RESULT="$?"
		if [ "$RESULT" != "0" ]; then
			exit 0
		fi
		"""

