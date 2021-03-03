import sys
Polymorphism_file = sys.argv[1]

#CCA     NZ_CP048438.1   2299440 2299440 NZ_AP021898.1   2266820 2266820 NZ_RDBQ01000001.1       209270  209270  Polymorphism
#NZ_CP048438.1   2751716 Akkermansia_muciniphila_strain_JCM_30893_chromosome,_complete_genome    G       A
dic = {}
with open(Polymorphism_file)  as F:
	for line in F:
		if "SNP" in line: continue
		l = line.rstrip().split()
		Variants = l[0]
		Contigs = len(l[0])
		n = 1
		for Entry in range(Contigs):
			Ref = Variants[Entry]
			Chr = l[n]
			Pos = l[n+1]
			for I in Variants:
				if I != Variants[Entry]: 
					Alt = I
					if Chr not in dic: dic[Chr] = {}
					if Pos not in dic[Chr]: dic[Chr][Pos] = {}
					if Alt not in dic[Chr][Pos]: dic[Chr][Pos][Alt] = ""
					dic[Chr][Pos][Alt] = Ref
			n+=3
output = sys.argv[2]
with open(output, "w") as F:
	for Chr in dic:
		for Pos in dic[Chr]:
			for Alt in dic[Chr][Pos]:
				Entry = "{Chr}\t{Pos}\t{Ref}\t{Alt}\n".format(Chr= Chr, Pos= Pos, Ref=dic[Chr][Pos][Alt], Alt= Alt)
				F.write(Entry)
