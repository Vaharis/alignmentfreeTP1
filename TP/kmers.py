
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
    encodage = {"A": 0, "C": 1, "T": 2, "G": 3}

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

            yield min(kmer1, kmer2)
