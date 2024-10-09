# TP1 - PHYG : Comparaison de séquences sans alignement


### Students : Océane LI et Adam BOUMESSAOUD


Le but de ce TP est de comparer 5 espèces bactériennes entre-elles, en se basant sur leur distance de Jaccard - paire par paire - qu'on va calculer.


## Matrice des distances (Jaccard)

Voici la matrice de Jaccard que nous avons obtenue :

| Bactéries          | GCA_000069965.1 | GCA_000008865.2 | GCA_030271835.1 | GCA_000013265.1 | GCA_000005845.2 |
|----------------------|------------------|------------------|------------------|------------------|------------------|
| GCA_000069965.1      | 1.0              | 0.002313787366673418 | 0.03113213059811918 | 0.002437015149030606 | 0.0025674665312398607 |
| GCA_000008865.2      | 0.002313787366673418 | 1.0              | 0.0023179077906273736 | 0.3070497016021145  | 0.43648444614279863  |
| GCA_030271835.1      | 0.03113213059811918 | 0.0023179077906273736 | 1.0              | 0.002433885647426537 | 0.0025764685291802055 |
| GCA_000013265.1      | 0.002437015149030606 | 0.3070497016021145  | 0.002433885647426537 | 1.0              | 0.3410085892813939  |
| GCA_000005845.2      | 0.0025674665312398607 | 0.43648444614279863  | 0.0025764685291802055 | 0.3410085892813939  | 1.0              |


Le temps de calcul pour avoir cette matrice est d'environ 50 secondes (vérifiée à l'aide de la fonction `time` mais peut changer en fonction de l'appareil).

## Analyse de la matrice Jaccard

La matrice ci-dessus présente les indices de Jaccard entre les séquences d'ADN de différentes bactéries. Chaque valeur représente la similarité entre les k-mers des bactéries une à une. Une valeur proche de 1.0 indique une grande similarité, tandis qu'une valeur proche de 0.0 indique une faible similarité.

- Espèces les plus différentes : **GCA_000069965.1** et **GCA_000008865.2**, avec un score de 0.002313787366673418 (soit seulement 0.23% de similarité).

- Espèces les plus similaires, **GCA_000008865.2** et **GCA_000005845.2**, avec un score de 0.43648444614279863 (soit ~43.65% de similarité).



## Méthodes implémentées


Nous avons construit une fonction `stream_kmers` qui génère des k-mers à partir des séquences d'ADN. Elle encode les nucléotides en valeurs entières et produit les k-mers et leurs compléments. Les k-mers sont représentés sous forme d'entiers pour accélérer le traitement.

Ensuite, la fonction `create_index` construit un dictionnaire où chaque clé est un k-mer et chaque valeur est le nombre d'occurrences de ce k-mer. Nous l'utilisons pour créer l'index de l'une des deux séquences à analyser, puis nous étudions toutes les paires possibles contenant cette séquence. Cela permet de compter efficacement les k-mers lors du calcul de l'indice de Jaccard.

Enfin, la fonction `jaccard` calcule l'indice de Jaccard entre les k-mers d'une première séquence (fichier A) et d'une seconde séquence (fichier B). Elle prend directement l'index du fichier A en argument. Elle maintient un compteur d'intersection et une somme de l'union des k-mers pour calculer la similarité.


**Conclusion :** Ce TP fournit un aperçu des relations entre divers ADN de bactéries à travers l'analyse des k-mers. La matrice de Jaccard permet d'identifier les espèces similaires et pourrait servir de base pour des analyses plus approfondies sur les relations phylogénétiques et évolutives.