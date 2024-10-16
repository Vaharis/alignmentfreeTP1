# TP2 - PHYG : Optimisation de la comparaison de séquences sans alignement


### Students : Océane LI et Adam BOUMESSAOUD


Le but de ce TP est de comparer 5 espèces bactériennes entre-elles, en se basant sur leur distance de Jaccard - paire par paire - qu'on va calculer. Afin d'optimiser le temps de calcul, nous avons utilisé une technique de hachage, pour à la fin ne comparer que des sketchs de taille *s* = 1000 kmers ici, plutôt que de comparer la totalité des kmers entre eux. Les matrices obtenues sont disponibles en format *.csv*.


## Matrice des distances de Jaccard (Bactéries)

Voici la matrice de Jaccard que nous avons obtenue :

| Bactéries       |   GCA_000005845.2 |   GCA_000008865.2 |   GCA_000013265.1 |   GCA_000069965.1 |   GCA_030271835.1 |
|-----------------|-------------------|-------------------|-------------------|-------------------|-------------------|
| GCA_000005845.2 |        1          |        0.37457    |        0.302932   |        0.00452034 |        0.00452034 |
| GCA_000008865.2 |        0.37457    |        1          |        0.276324   |        0.00351229 |        0.00351229 |
| GCA_000013265.1 |        0.302932   |        0.276324   |        1          |        0.00401606 |        0.00401606 |
| GCA_000069965.1 |        0.00452034 |        0.00351229 |        0.00401606 |        1          |        0.0465725  |
| GCA_030271835.1 |        0.00452034 |        0.00351229 |        0.00401606 |        0.0465725  |        1          |

Temps de calcul: ~ 20sec

## Matrice des distances de Jaccard (Mammifères)

| Species          |   GCA_029289425.3 |   GCF_000001405.40 |   GCF_000001635.27 |
|------------------|-------------------|--------------------|--------------------|
| GCA_029289425.3  |        1          |         0.473839   |         0.00401606 |
| GCF_000001405.40 |        0.473839   |         1          |         0.00401606 |
| GCF_000001635.27 |        0.00401606 |         0.00401606 |         1          |

Temps de calcul: ~1h15

## Analyse de la matrice Jaccard

Les matrices présentent les indices de Jaccard entre les séquences d'ADN de différentes espèces. Chaque valeur représente la similarité entre les k-mers des espèces une à une. Une valeur proche de 1.0 indique une grande similarité, tandis qu'une valeur proche de 0.0 indique une faible similarité.
En ce qui concerne les bactéries, on obtient des résultats similaires à l'implémentation précédente: les 3 espèces les plus proches sont toujours GCA_000005845.2, GCA_000008865.2 et GCA_000013265.1.
Dans la matrice des mammifères, on observe une similarité relativement haute entre l'Homo Sapiens (GCF_000001405.40) et le Pan paniscus (GCA_029289425.3) avec un indice de 0.47 environ. Le Mus Musculus cependant est relativement éloigné des deux autres espèces.
La matrice comparant les bactéries aux mammifères est disponible en format *.csv*, il n'y a aucune similarité visible avec les indices de Jaccard obtenus.

## Méthodes implémentées

Nous avons utilisé l'algorithme de `bottom_minhach` dans le but de ne sélectionner que les *s* plus petits kmers (les kmers sont codés sous forme d'entiers). En faisant cela, nous avons introduit un biai lié à l'ordre lexical des nucléotides: A=0, C=1, T=2, G=3. 

Afin d'éviter la sélection biaisée de kmers commençant par de nombreuses Adénines, nous avons utilisé une fonction de hachage `xorshift` mélangeant ainsi les kmers de manière prévisible (deux kmers identiques passant par cette fonction restent identiques en sortant de celle-ci).

Nous avons modifié la fonction `jaccard`, elle permet maintenant de calculer l'indice de Jaccard à partir de deux listes de nombres (nombres correspondants aux *s* kmers choisis)

Nous avons retiré des séquences étudiées les régions soft-masked (nucléotides écrites en minuscule) car elles correspondent à des régions répétitives du génome. Nous nous sommes également affranchis des nucléotides ambigus car ces derniers ne correspondent qu'à une petite partie du génome et que procéder ainsi permet de grandement améliorer la vitesse du programme.

**Conclusion :** Lors de ce TP nous avons implémenté une méthode heuristique permettant de calculer l'indice de jaccard entre espèces de façon beaucoup plus rapide que lors du TP1. Nous avons eu l'occasion tester notre travail sur des séquences de taille conséquente (3Go et +), ce qui nous a permis de bien comprendre l'importance de l'optimisation de notre code.
