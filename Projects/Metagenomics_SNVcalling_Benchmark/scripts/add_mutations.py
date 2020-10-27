import sys
from Bio import SeqIO
import random
from subprocess import call
from pathlib import Path
import argparse

random.seed(2)

parser = argparse.ArgumentParser()
parser.add_argument('--fasta', '-f', action = 'store', default = '', dest = 'FASTA',required=True)
parser.add_argument('--N_mutations', '-m', action = 'store', default = '', dest = 'N_mutations', required=True)
parser.add_argument('--output', '-o', action = 'store', default = '', dest = 'OUT', required=True)
parser.add_argument('--output_m', '-om', action = 'store', default = '', dest = 'OUTm', required=True)


FASTA = parser.parse_args().FASTA
N_mutations = parser.parse_args().N_mutations
output = parser.parse_args().OUT
output_m = parser.parse_args().OUTm
if "%" in N_mutations:
	Proportion = True
	N_mutations_p = float(N_mutations.strip("%"))
else:
	Proportion = False
	N_mutations = int(N_mutations)

mutations_options = {"A":["T","C","G"], "T":["A","C","G"], "C":["T","A","G"], "G":["T","C","A"], "N":["A","T","G","C"]}
dic_mutations = {}
with open(output, "w") as OUT: pass
with open(FASTA, "r") as File:
	for record in SeqIO.parse(File, "fasta"):
		NAME = record.id
		SEQ =  record.seq.tomutable()
		Length_SEQ = len(SEQ)
		complete_name =  record.description
		chrom = "_".join(complete_name.split(" ")[1:])

		dic_mutations[NAME] = {}
		if Proportion == True: N_mutations = int(Length_SEQ * N_mutations_p/100)
		Position_mutated = random.choices(range(0, Length_SEQ-1), k=N_mutations)
		for Position in Position_mutated:
			Prev = SEQ[Position]
			if Prev not in ["A","T","G","C"]: New = random.choice(["A","T","G","C"])
			else: New = random.choice(mutations_options[Prev])
			dic_mutations[NAME][Position] = [Prev, New, chrom]
			SEQ[Position] = New
		with open(output, "a") as OUT:
			OUT.write(">"+complete_name+"\n"+str(SEQ)+"\n")

with open(output_m,"w") as O:
	for Genome in dic_mutations:
		for mutation in dic_mutations[Genome]:
			Line = "\t".join([Genome, str(mutation+1),dic_mutations[Genome][mutation][2], dic_mutations[Genome][mutation][0], dic_mutations[Genome][mutation][1]]) + "\n"
			O.write(Line)			

