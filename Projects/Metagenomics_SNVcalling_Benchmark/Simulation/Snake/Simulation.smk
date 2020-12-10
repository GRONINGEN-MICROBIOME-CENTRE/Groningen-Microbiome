def Polymor(Working_dir = Working_dir, design = design, IDs = dic_ref):
	To_add = expand("{working_dir}Simulation/Variants/{ID}_mutations.tsv",working_dir=Working_dir, ID= IDs)
	if design != "uniref":
		To_add2 = expand("{working_dir}Simulation/Variants/Polymorphisms/{ID}.snp",working_dir=Working_dir, ID= IDs)
		To_add.extend(To_add2)
	return(To_add)
def Get_Ref(ID):
	return(dic_ref[ID][0])
 


Mauve="/groups/umcg-lld/tmp04/umcg-sandreusanchez/Pipelines/SNV_calling/Simulation/Alignments/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve"


rule Find_polymorphism:
	input: Ref = lambda wildcards: Get_Ref(wildcards.ID),
	output: 
		Alignment = temp(Working_dir + "Simulation/Variants/Polymorphisms/{ID}.fa"),
		SNP = Working_dir + "Simulation/Variants/Polymorphisms/{ID}.snp" 
	shell:
		"""
		ml Java/11.0.2
		MAUVE_DIR="/groups/umcg-lld/tmp04/umcg-sandreusanchez/Pipelines/SNV_calling/Simulation/Alignments/mauve_snapshot_2015-02-13/"
		find "$MAUVE_DIR" -name \\*.jar -print0 >> CHECKME
		export CLASSPATH=""
		
		export CLASSPATH="$(find "$MAUVE_DIR" -name \\*.jar -print0 | tr '\\0' :)$CLASSPATH" 
		
        	{Mauve} --output {output.Alignment} {input.Ref}
		java org.gel.mauve.analysis.SnpExporter -f {output.Alignment}  -o {output.SNP}	
		"""
rule Introduce_mutations:
	input: Ref = lambda wildcards: Get_Ref(wildcards.ID),
	output:
		Fasta = Working_dir + "Simulation/Reference/{ID}.fa",
		Mutations = Working_dir + "Simulation/Variants/{ID}_mutations.tsv"
	params:
		mutations = config["Number_mutations"],
		py = "/home/umcg-sandreusanchez/.conda/envs/myenv/bin/python3.6" #Python 3.6 required for random.choices ; biopython needed
	shell:
		"if {design}==uniref; then ;\n"
		"{params.py} scripts/add_mutations.py --fasta {input.Ref} --N_mutations {params.mutations} --output {output.Fasta} --output_m {output.Mutations} --N_strains 1 ;\n"
		"else;\n"
		"{params.py} scripts/add_mutations.py --fasta {input.Ref} --N_mutations {params.mutations} --output {output.Fasta} --output_m {output.Mutations} --N_strains 2;\n"
		"fi"

rule Simulate:
	input: expand(Working_dir+"Simulation/Reference/{ID}.fa",ID=dic_ref)
	output: 
		#Ref = temp(Working_dir + "Simulation/Reference/Reference.fa"),
		flag = touch(Working_dir +"Simulations_completed.flag"),
		Read_f =  touch(Working_dir +"Simulation/FQ/simulated_R1.fastq"),
		Read_r = touch(Working_dir +"Simulation/FQ/simulated_R2.fastq"),
	params:
		n_reads = config["Read_number"],
		PrefixOut = temp(Working_dir + "Simulation/FQ/simulated"),
		Abundance = config["abundance_file"]
	shell:
		"""
		module load  Miniconda3/4.4.10 ; source activate  /groups/umcg-lld/tmp04/umcg-sandreusanchez/Conda_env/InSilicoSeq

		if [ ! -d "{Working_dir}Simulation/FQ/" ]; then  mkdir {Working_dir}Simulation/FQ/ ; fi
		if [ {params.Abundance} != "NA" ]
		then
			set +e
			iss generate --n_reads {params.n_reads} --draft {input} --model hiseq --output {params.PrefixOut} --seed 5 --abundance_file {params.Abundance} --cpus 1
			exitcode=$?
			if [ $exitcode -eq 1 ]
			then
			exit 0
			fi
		else 
			iss generate --n_reads {params.n_reads} --draft {input} --model hiseq --output {params.PrefixOut} --seed 5 --cpus 1
		fi
		"""
		


rule Combine_polymorphisms:
	input: Polymor()
	output:
		Mutations = Working_dir +"Simulation/Mutations_introduced.tsv",
		Polymorphisms = Working_dir +"Simulation/Polymorphisms.tsv",
		Polymorphisms2 = Working_dir +"Simulation/Polymorphisms_formatted.tsv"
	params:
		design = config["design"]
	shell:
		"""
		cat {Working_dir}Simulation/Variants/*_mutations.tsv > {output.Mutations}
		if [ {params.design} == "uniref" ]
		then
			touch {output.Polymorphisms2} ; touch {output.Polymorphisms}
		else
			cat {Working_dir}Simulation/Variants/Polymorphisms/*.snp > {output.Polymorphisms}
			python scripts/Format_polymorphism_mauve.py {output.Polymorphisms} {output.Polymorphisms2}
		fi
		"""
