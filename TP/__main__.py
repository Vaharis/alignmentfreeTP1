from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str
import numpy as np
import pandas as pd
import time
start_time = time.time()

def create_index(fileA, k):
	"""
	Prépare l'index contenant en clé un kmer et en valeur son nombre d'occurences à partir d'une séquence et d'une longueur k.
	----------
	Params:
		- fileA: séquences ADN du fichier 1.
		- k: taille des k-mers.
	----------
	Returns:
		- index: Le dictionnaire avec en clé un kmer et en valeur son nombre d'occurences.
	"""
	# Creating the index based on the first sequence
	index = {}
	for kmer in stream_kmers(fileA, k):
		if kmer not in index:
			index[kmer] = 1
		else:
			index[kmer] +=1
	return index

def jaccard(indexA, fileB, k):
	"""
	Calcule l'indice de Jaccard entre deux fichiers de séquences d'ADN.
	----------
	Params:
		- indexA: indexe des séquences ADN du fichier 1.
		- fileB: séquences ADN du fichier 2.
		- k: taille des k-mers.
	----------
	Returns:
		- jaccard_index: Indice de Jaccard entre les 2 fichiers.
	"""

    # Calculation of the jaccard distance
	union = sum(indexA.values())
	intersection = 0
	
	for kmer in stream_kmers(fileB, k):
		if kmer in indexA and indexA[kmer] > 0:
			intersection += 1
			indexA[kmer] -= 1
		else:
			union += 1
	
	return intersection/union


if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    jaccard_matrix = np.zeros((len(filenames), len(filenames)))
    np.fill_diagonal(jaccard_matrix, 1)

    for i in range(len(files)):

		# On créé l'index de la séquence A ici car il sera le même pour tous les B suivantes.
        indexA = create_index(files[filenames[i]], k)

        for j in range(i + 1, len(files)):
			
            jaccard_score = jaccard(indexA.copy(), files[filenames[j]], k) # Puisque l'index sera modifié dans la fonction, on doit créer une copie avec index.copy()
            jaccard_matrix[i][j] = jaccard_score
            jaccard_matrix[j][i] = jaccard_score

            print(filenames[i], filenames[j], jaccard_score)

    print("\n",pd.DataFrame(data=jaccard_matrix, columns = filenames, index=filenames))
    print("\n--- %s seconds ---" % (time.time() - start_time))
