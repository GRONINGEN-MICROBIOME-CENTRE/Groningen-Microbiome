rule Coverage:
	input:
		BAM=  Working_dir+"Alignment/Aligned/{ID}_clean.bam",
		REF= lambda wildcards: Get_Ref(wildcards.ID)
	output:
		cov = Working_dir+"Ref_variants/Coverage/{ID}.cov"
	shell:
		"""
		ml SAMtools/1.9-foss-2015b BCFtools/1.7-foss-2015b
		samtools mpileup --fasta-ref {input.REF} -A {input.BAM} > {output.cov}
		"""

rule Find_Covered_variants:
	input:
		Polymor = Working_dir +"Simulation/Polymorphisms2.tsv",
		Mutation =  Working_dir +"Simulation/Mutations_introduced.tsv",
		Covered_positions = Working_dir+"Ref_variants/Coverage/{ID}.cov"
	output:
		cov = Working_dir+"Ref_variants/{ID}_covMutations",
		mutations_c = temp(Working_dir+"tmp/mutations_covered_{ID}"),
		poly_c =  temp(Working_dir+"tmp/polym_covered_{ID}"),
	shell:
		"""
		bash scripts/Check_mutations.sh {input.Polymor} {input.Mutation} {output.mutations_c} {output.poly_c} {output.cov} {input.Covered_positions}
		"""
rule Check_if_called:
	input:
		Covered =  Working_dir+"Ref_variants/{ID}_covMutations",
		Called =  Working_dir+"Calls/{tool}/{ID}.tsv",
		coverage =   Working_dir+"Ref_variants/Coverage/{ID}.cov"
	
	output:
		tempcall = temp(Working_dir+"Calls/{tool}/{ID}_unique.tsv"),
		Tool_out = Working_dir+"Statistics_calling/{tool}.{ID}.tsv"
	shell:
		"""
		sort {input.Called} | uniq > {output.tempcall}	
		ml Biopython/1.65-foss-2015b-Python-3.4.1 ; python scripts/Classify_Variants.py --ref {input.Covered} --calls {output.tempcall} --covered {input.coverage} --output {output.Tool_out} --threshold 0

		"""

rule ROC_MultipleThresholds:
	input:
		Covered =  Working_dir+"Ref_variants/{ID}_covMutations",
		Called =  Working_dir+"Calls/{tool}/{ID}.tsv",
		coverage =   Working_dir+"Ref_variants/Coverage/{ID}.cov"
	output:
		tempcall = temp(Working_dir+"Calls/{tool}/{ID}_unique2_{Threshold}.tsv"),
		temporal_output1 =temp( Working_dir + "tmp/{tool}.{ID}.{Threshold}"),
		temporal_output2 = temp(Working_dir + "tmp/{tool}.{ID}.{Threshold}_stats")
	shell:
		"sort {input.Called} | uniq > {output.tempcall} ;\n"
		"awk '$5 > {wildcards.Threshold}' {output.tempcall} > {output.temporal_output1} ;\n"
		"ml Biopython/1.65-foss-2015b-Python-3.4.1 ; python scripts/Classify_Variants.py --ref {input.Covered} --calls {output.temporal_output1} --covered {input.coverage} --output {output.temporal_output2} --threshold {wildcards.Threshold}"

def ROC_thresholds(wildcards, DIR= Working_dir):
	files = []
	for threshold in range(0, 5000, 20):
		I = "{DIR}tmp/{tool}.{ID}.{Threshold}_stats".format(DIR=DIR, tool = wildcards.tool, ID=wildcards.ID, Threshold=threshold)
		files.append(I)
	return(files)
rule ROC:
	input:
		lambda wildcards: ROC_thresholds(wildcards)
	output:
		ROC = Working_dir+"Statistics_calling/{tool}.{ID}.tsvROC"
	params:
		tempfolder = Working_dir + "tmp/{tool}_{ID}"
	run:    
		for F in Inputs:
			#TP      TN      FP      FN      Taxa    Tool    Length_genome   Benchmark_Design
			with open(F) as FILE:
				Thresh = F.split("_")[-2]
				for line in FILE:
					if "TP" in  line: continue
					l = line.rstrip().split()
					TPR = float(l[0])/(float(l[0])+float(l[2]))
					FPR = float(l[2])/(float(l[0])+float(l[2]))
					Taxa = l[4]
					with open(output.ROC, "a") as O:
						O.write("\t".join([str(TPR), str(FPR), str(Thresh), Taxa])+"\n")
