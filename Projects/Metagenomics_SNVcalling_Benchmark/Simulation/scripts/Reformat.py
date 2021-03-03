import sys
from pathlib import Path
from subprocess import call
import argparse
#Script meant to generate a single Output file
#Inputs: 
#	vcf files, program-specific format


parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', action = 'store', dest = 'input')
parser.add_argument('--program', '-p', action = 'store', dest = 'program')
parser.add_argument('--output', '-o', action = 'store', dest = 'Out') #, choices = set(('Unistrain','Multistrain','Real')))
args = parser.parse_args()

File = args.input
Program = args.program
Output =args.Out

def Add_to_file(Info, Output=Output):
	Output = Path(Output)
	if not Output.parent.exists():
		call("mkdir {P}".format(P=Output.parent), shell=True)
	with open(str(Output),"a") as F:
		if Info[-1] == ".": Info = Info[0:-1]
		F.write("\t".join(Info) + "\n")


def Format_vcf_old(File, Output):
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER
	with open(File) as F:
		for line in F:
			if line[0] == "#": continue
			l = line.rstrip().split()
			if len(l) < 5: continue
			Chrom = l[0]
			Pos = l[1]
			Ref = l[3]
			Alt = l[4]
			if l[5] != ".": QUAL = l[5]
			else: QUAL = None
			if "," in l[4]:
				for Alt in l[4].split(","):
					if not QUAL is None:  Add_to_file([Chrom, Pos, Ref, Alt, QUAL], Output)
			else:
				Add_to_file([Chrom, Pos, Ref, Alt], Output)
def Format_vcf(File, Output):
	#Bacteroides_uniformis   222     .       TGCTG   CGCTA,GGCTG,TGATG
	with open(File) as F:
		for line in F:
			if line[0] == "#": continue
			l = line.rstrip().split()
			if len(l) < 5: continue
			Chrom = l[0] ; Pos = l[1] ; Ref = l[3] ; Alt = l[4] ; QUAL = l[5]
			if "," in l[4]:
				for Alt in l[4].split(","):
					if len(Alt) > 1:
						for n in range(len(Alt)):
							if n >= len(Ref): continue
							if Alt[n] != Ref[n]:
								Add_to_file([Chrom, str(int(Pos)+n), Ref[n], Alt[n], QUAL], Output)
					else: Add_to_file([Chrom, Pos, Ref, Alt, QUAL], Output)
			else:
				if len(Alt) > 1:
					for n in range(len(Alt)):
						if n >= len(Ref): continue
						if Alt[n] != Ref[n]:
							Add_to_file([Chrom, str(int(Pos)+n), Ref[n], Alt[n], QUAL], Output)
				else:
					if len(Ref) > len(Alt):
						if Ref[0] == Alt: continue
						Add_to_file([Chrom, Pos, Ref[0], Alt,QUAL], Output)
					else:
						Add_to_file([Chrom, Pos, Ref, Alt,QUAL], Output)
					
def Format_metasnv(File, Output):
	#NZ_CP048438.1   -       478857  A       2|2|0|0 4|T|.|2|2|0|0
	#Faecalibacterium_prausnitzii    -       259127  T       4       4|A|.|4,
	with open(File) as F:
		for line in F:
			l = line.rstrip().split()
			Chrome=l[0] ; Pos=l[2] ; Ref=l[3]
			if "," in l[5]:
				for parts in l[5].split(","):
					Variant = parts.split("|")[1]
					To_write = [Chrome, Pos, Ref, Variant]
					Add_to_file(To_write, Output)
			else:
				Variant = l[5].split("|")[1]
				To_write = [Chrome, Pos, Ref, Variant]
				Add_to_file(To_write, Output)

def Format_instrain(File, Output):
	#scaffold	position	refBase	A	C	T	G	conBase	varBase	allele_count	cryptic	baseCoverage	varFreq	conFreq	refFreq
	with open(File) as F:
		for line in F:
			l = line.rstrip().split()
			if len(l) <= 4: continue
			Chr = l[0]
			Pos = l[1]
			Ref = l[2]
			if "scaffold" in line: continue
			if l[3] != "0":
				Alt = "A"
				if Alt != Ref: 
					Add_to_file([Chr, Pos, Ref, Alt], Output)
			if l[4] != "0":
				Alt = "C"
				if Alt != Ref:
					Add_to_file([Chr, Pos, Ref, Alt], Output)
			if l[5] != "0":
				Alt = "T"
				if Alt != Ref:
					Add_to_file([Chr, Pos, Ref, Alt], Output)
			if l[6] != "0":
				Alt = "G"
				if Alt != Ref:
					Add_to_file([Chr, Pos, Ref, Alt], Output)



with open(Output, "w"): pass


if Program=="metasnv": #/groups/umcg-lld/tmp04/umcg-lchen/metaSNV/result_Akkermansia_muciniphilancbi/snpCaller/called_SNPs
	Format_metasnv(File,Output)
elif Program == "instrain":
	Format_instrain(File,Output)
else:
	Format_vcf(File, Output)
