import random

def kmer2str(val, k):
    """ Transform a kmer integer into its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)

def xorshift(val):
	"""
	Fonction de hachage.
	----------
	Params:
		- val: int;  la valeur à hacher
	----------
	Returns:
		- val: int; la valeur hachée
	"""
	val ^= val << 13
	val &= 0xFFFFFFFFFFFFFFFF
	val ^= val >> 7
	val ^= val << 17
	val &= 0xFFFFFFFFFFFFFFFF
	return val
    
def stream_kmers(text, k):
    """
    Generates k-mers and their canonical form from a list of sequence texts.
    ----------
    Params:
        - text (list[str]): chaque str est une séq d'ADN.
        - k (int): taille des k-mers à générer.
    ----------
    Returns:
        - le k-mer canonique à chaque itération.
    """

    # Masque binaire pour garder uniquement les bits correspondant à un k-mer de taille k
    mask_bin = (1 << (2 * k)) - 1
	
    # encodage : dict{} ; mappe chaque nucléotique à un int
    encodage = {
                 'A': 0, 'C': 1, 'T': 2, 'G': 3, 'R': 3, 'Y': 1, 'K': 2, 'M': 0, 
                 'S': 1, 'W': 2, 'B': 3, 'D': 0, 'H': 1, 'V': 3, 
                 'a': 0, 'c': 1, 't': 2, 'g': 3, 'r': 3, 'y': 1, 'k': 2, 'm': 0, 
                 's': 1, 'w': 2, 'b': 3, 'd': 0, 'h': 1, 'v': 3
                }
    encodage = {'A': 0, 'C': 1, 'T': 2, 'G': 3}
    for seq in text:  # on itère sur chaque seq
        
        # kmer1, kmer2 ; représentation binaire du k-mer et de son complément inverse
        kmer1, kmer2 = 0, 0
        
        # on initialise les k-mers avec les premiers k-1 nucléotides
        for letter in seq[: k - 1]:
            nucl = encodage[letter]  # int value du nt
            comp = (nucl + 2) & 3  # calcul du nt complémentaire

            # à chaque ajout de nv nt, les kmers sont décalés de 2 bits à gauche
            kmer1 = ((kmer1 << 2) + nucl) & mask_bin  # nt
            kmer2 = ((kmer2 << 2) + comp) & mask_bin  # complément
        
        
        # Reste de la séquence :
        for letter in seq[k - 1 :]:
            nucl = encodage[letter]
            comp = (nucl + 2) & 3  # nt complément
            comp <<= 2 * k - 2  # ajustement pour le complément

            kmer1 = ((kmer1 << 2) + nucl) & mask_bin
            kmer2 = ((kmer2 >> 2) + comp) & mask_bin

            yield min(xorshift(kmer1), xorshift(kmer2))
