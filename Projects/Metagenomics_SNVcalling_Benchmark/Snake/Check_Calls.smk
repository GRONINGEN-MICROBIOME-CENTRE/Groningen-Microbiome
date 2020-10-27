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
rule ROC:
        input:
                Covered =  Working_dir+"Ref_variants/{ID}_covMutations",
                Called =  Working_dir+"Calls/{tool}/{ID}.tsv",
                coverage =   Working_dir+"Ref_variants/Coverage/{ID}.cov"    
        output:
                tempcall = temp(Working_dir+"Calls/{tool}/{ID}_unique2.tsv"),
                ROC = Working_dir+"Statistics_calling/{tool}.{ID}.tsvROC"
        params:
                tempfolder = Working_dir + "tmp/{tool}_{ID}"
        run:    
                shell("sort {input.Called} | uniq > {output.tempcall}")
                list_files =[] 
                for Threshold in range(0,305, 5): 
                        temporal_output1 = params.tempfolder + "_"+ str(Threshold)
                        temporal_output2 = params.tempfolder + "_"+ str(Threshold) + "_stats"
                        shell("awk '$5 > {Threshold}' {output.tempcall} > {temporal_output1}")   
                        shell("ml Biopython/1.65-foss-2015b-Python-3.4.1 ; python scripts/Classify_Variants.py --ref {input.Covered} --calls {temporal_output1} --covered {input.coverage} --output {temporal_output2} --threshold {Threshold}")
                        shell("rm {temporal_output1}")
                        list_files.append(temporal_output2)
                for F in list_files:
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
