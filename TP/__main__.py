from TP.loading import load_directory # type: ignore
from TP.kmers import stream_kmers, kmer2str, xorshift # type: ignore
import numpy as np
import pandas as pd
import time
import heapq
start_time = time.time()

def bottom_minhash(seq, k, s):
	"""
	Crée un sketch (une liste) contenant les s plus petits kmers hachés canoniques de taille k d'une séquence.
	Les kmers sont ordonnés de manière croissante (ils sont codés sous forme d'entiers après le hachage). 
	----------
	Params:
		- seq: list;  la séquence à analyser
		- k: int; taille des k-mers.
		- s: int; taille du sketch
	----------
	Returns:
		- sorted_sketch: la liste à retourner
	"""

	# Création du sketch
	sketch = [-np.inf for _ in range(s)]
	heapq.heapify(sketch)

	for kmer in stream_kmers(seq,k):
		if -kmer > sketch[0]:
			heapq.heappushpop(sketch, -kmer)

	# Tri du sketch
	sorted_sketch = []
	while sketch:
		sorted_sketch.append(heapq.heappop(sketch))
	
	return sorted_sketch

def jaccard(sketchA, sketchB):
	"""
	Calcule l'indice de Jaccard entre deux listes.
	----------
	Params:
		- sketchA: La liste 1.
		- sketchB: La liste 2.
	----------
	Returns:
		- intersection/union : float correspondant à l'indice de jaccard
	"""
	# i, j: int; correspondent aux indices d'itération sur la 1ère et 2e liste respectivement
	i = 0
	j = 0

	union = 0
	intersection = 0

	while i < len(sketchA) and j < len(sketchB):

		if sketchA[i] < sketchB[j]:
			union += 1
			i += 1

		elif sketchA[i] > sketchB[j]:
			union += 1
			j += 1

		elif sketchA[i] == sketchB[j]:
			union += 1
			intersection += 1
			i += 1
			j += 1

	union += abs(len(sketchA) - i + len(sketchB) - j)
	return intersection/union


if __name__ == "__main__":
	print("Computation of Jaccard similarity between files")
	

	# Lecture des séquences
	files = load_directory("data")
	filenames = list(files.keys())

	# Initialisation des paramètres
	k = 21		# Taille des kmers
	s = 1000	# Tailel des sketchs

	# Création des sketchs.
	print("\nCreation of the sketches.")
	sketch_list = []
	for i in range(len(filenames)):
		sketch_list.append(bottom_minhash(files[filenames[i]], k, s))
		print(f"{i+1}/{len(filenames)} sketches done")
	print()

	# Calcul des indices de jaccard
	jaccard_matrix = np.zeros((len(filenames), len(filenames)))		# Matrice contenant les résultats
	np.fill_diagonal(jaccard_matrix, 1)								# L'indice de jaccard de deux séquences identiques est de 1

	for i in range(len(filenames)):
		for j in range(i+1, len(filenames)):
				jaccard_score = jaccard(sketch_list[i], sketch_list[j])

				jaccard_matrix[i][j] = jaccard_score
				jaccard_matrix[j][i] = jaccard_score

				print(filenames[i], filenames[j], jaccard_score)

	jaccard_df = pd.DataFrame(data=jaccard_matrix, columns = filenames, index=filenames)
	#jaccard_df.to_csv("jaccard_df.csv", index_label = "Species")
	
	print("\njaccard_df.csv created")
	print("\n--- %s seconds ---" % (time.time() - start_time))
