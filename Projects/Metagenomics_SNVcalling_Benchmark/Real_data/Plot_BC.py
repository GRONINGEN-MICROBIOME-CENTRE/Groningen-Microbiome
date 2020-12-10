import numpy as np

Distance_matrix = "result_BC_distmat_86samples_318species_hmp.txt"
from scipy import spatial,stats,cluster
import matplotlib.pyplot as plt

Intra = []
Inter = []

line = 0 
Matrix = False
with open(Distance_matrix) as F:
	for l in F:
		if line == 0: 
			Header = l.rstrip().split()
			line +=1
			continue
		Profile = np.array(l.rstrip().split()[1:])
		if isinstance(Matrix, (bool)): Matrix= Profile
		else: Matrix = np.vstack((Matrix, Profile))

print(Matrix.shape)


Intra_samples = []
for N in range(len(Header)):
	H1 = Header[N].replace("_F","")
	for N2 in range(len(Header)):
		H2 = Header[N2].replace("_F","")
		if N != N2 and H1==H2:
			D = str(N) + "-" + str(N2)
			Intra_samples.append(D)
		
Distance_triangular = np.triu(Matrix)
i1= 0
Intra_individual = []
Inter_individual = []
for Row in Distance_triangular:
	i2 =0
	for column in Row:
		if column == "": i2 +=1 ; continue
		Position =  str(i1) +"-" + str(i2)
		if column == 0: i2 +=1 ; continue
		if Position in Intra_samples: 
			Intra_individual.append(float(column))
			i2 +=1
			continue
		if i1 == i2:
			i2 +=1 ; continue
		Inter_individual.append(float(column))			
		i2 +=1
	i1+=1


print("Performing wilcox test between inter and intra individual distanes")
w, p = stats.mannwhitneyu(Intra_individual, Inter_individual)
print("Pvalue=" + str(p))

M = abs(np.mean(np.array(Intra_individual)) - np.mean(np.array(Inter_individual)))

print("Mean difference=" + str(M))


print("Plotting distributions")
data = [np.array(Intra_individual), np.array(Inter_individual)]

fig, ax = plt.subplots()
ax.boxplot(data)
positions = (1, 2)
labels = ("Intra-individual", "Inter-individual")
plt.xticks(positions, labels)
plt.savefig("Plots/Inter-vs-Intra_abundance_BC.png")	
