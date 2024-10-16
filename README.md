
# Alignment free - TP 1

## Télécharger la base de code du TP

Durant nos TPs, nous allons utiliser Git pour récupérer le TP, versionner son code et rendre le TP. 

### Je n'ai jamais utilisé Git

Si tu n'as jamais utilisé Git, ce TP est le bon moment pour démarrer. En bioinformatique, tous les logiciels produits sont hébergés sur des dépôts en ligne Git.
Pour te faire la main avec git, je te laisse parcourir [ce tutoriel qui t'expliquera les bases](https://git.goffinet.org/02-les-bases-de-git).

Pour ce TP, nous utiliserons Github comme hébergeur. Si vous n'avez pas de compte, il vous faudra en créer un pour pouvoir rendre le TP.
Une fois votre compte créé vous allez devoir enregistrer votre clé ssh sur le serveur. Dans les paramètres de votre compte, allez dans la rubrique "SSH and GPG keys" et suivez le lien qui explique comment entrer votre clé publique.

Vous êtes maintenant prêt.e à télécharger les bases de votre TP.

### Télécharger les bases pour le TP

* Je me rends à l'adresse du TP : https://github.com/yoann-dufresne/alignmentfreeTP1

* Je click sur le bouton "fork" en haut de la page pour créer une copie du TP au sein de mon compte. Je renomme la copie comme bon me semble et je valide.

* Sur la page où j'arrive, je récupère l'adresse à laquelle télécharger le code. Pour cela je clique sur "code" puis "ssh" et je copie l'adresse qui apparait dans l'interface.

* Enfin, sur mon ordinateur, à l'emplacement où je souhaite télécharger le TP, j'utilise la commande suivante **en remplaçant l'adresse git par celle que je viens de copier**.

```bash
    # Télécharge le dépôt de code initial. Le --recursive est nécessaire pour télécharger les sous-dépôts liés.
    git clone --recursive git@github.com:ada-lovelace/alignmentfreeTP1.git
```


## Objectifs du TP

Pour ce TP il faudra comparer un ensemble de 5 bactéries pour retrouver les familles présentes. L'objectif final est de produire une matrice des indices de Jaccard entre tous les échantillons.
**Attention** : Ce TP est une base pour le prochain. Si je ne le finis pas, je ne pourrai pas faire le suivant !

### Récupérer les séquences bactériennes

Pour télécharger les données un fichier dataset.reg qui contient l'ensemble des accessions ENA à télécharger à été inclus dans le dépôt.
Le logiciel seqdd qui permet de télécharger des données de séquences depuis plusieurs sources, a été clonné comme sous-module lorsque j'ai fait le git clone récursif. Il est facilement utilisable pour récupérer les fichiers fasta.
Pour cela, depuis la racine de mon projet je dois effectuer la liste des commandes suivantes:

```bash
    # Installe localement l'outil seqdd depuis le répertoire nommé seqdd
    python3 -m pip install --user seqdd
    # Initialise le registre de données avec le fichier .reg du tp
    seqdd init -r dataset.reg
    # Télécharge les données depuis le European Nucleotide Archive (ENA)
    seqdd download -d data/
```

Si cette procédure ne fonctionne pas, je peux toujours ouvrir le fichier .reg comme un fichier texte puis me rendre manuellement sur le [site web de l'ENA](https://www.ebi.ac.uk/ena/browser/home) et y télécharger manuellement les jeux de données listés.

### Composer le TP

Le code téléchargé fonctionne déjà mais ne calcule rien.
Pour le lancer, je dois appeler la commande :

```bash
    # Déclenche l'exécution du fichier __main__.py dans le module TP
    python3 -m TP
```

Si tout se passe bien, le code devrait charger tous les échantillons présents dans le répertoire data/ puis lancer une distance de Jaccard entre toutes les pairs d'échantillon.
La fonction Jaccard n'étant pas codée, toutes les valeurs affichées devraient être à 0.

Mon travail est de compléter les fonctions pour que les indices de Jaccard soient bons.


### Sauvegarder mon TP

En tant que bon.ne développeur.euse, mon TP est modifié et testé de manière incrémentale. Je m'assure à tout moment que ce que je viens de coder est correct.
Dès que j'ai ajouté une nouvelle fonctionnalité et que tout est ok, je la sauvegarde sur github de la manière suivante :

```bash
    # Ajout des fichiers modifiés à la liste des fichiers suivis par git
    git add TP/monfichier.py TP/mondeuxiemefichier.py
    # Packaging des fichiers modifiés et message explicatif de la nouvelle fonctionnalité
    git commit -m "Ajout de la fonctionnalité X, désormais le TP fait Y"
    # Envoie du package sur les serveurs github dans la branche main
    git push origin main
```

### Rendre mon TP

Le rendu attendu est un email à l'enseignant en indiquant l'adresse de mon dépôt git contenant le TP.
Je dois créer une branche "rendu_1" que je ne toucherai plus après avoir envoyé le mail. Voici la liste de commandes pour créer une branche à partir de main et l'envoyer sur le serveur :

```bash
    # Checkout permet d'aller sur la branche indiquée (ici ma_nouvelle_branche).
    # Le -b dit à git que la branche n'existe pas encore et qu'il faut la créer à partir du code actuel.
    git checkout -b ma_nouvelle_branche
    # Envoyer ma branche sur le serveur
    git push origin ma_nouvelle_branche
    # Revenir à ma branche principale
    git checkout main
```

Dans mon dépôt est attendu :
* Le code complété pour que les indices de Jaccard soient calculés
* Un fichier RENDU.md qui contient un rapport succinct sur mon rendu. Ce rapport doit contenir la matrice des distances entre échantillons, un commentaire sur ce que cette matrice m'apprend et une description des méthodes implémentées.

Le TP peut être rendu par bionôme. Je dois l'indiquer dans RENDU.md en précisant les prénoms et noms.


### Attendus pour le TP 1

Le code que vous avez produit lors de la première séance de TP doit s'exécuter en quelques secondes au maximum.
La matrice des identitées de Jaccard doit être très proche de:
|                 | GCA_000013265.1 | GCA_000005845.2 | GCA_000069965.1 | GCA_000008865.2 | GCA_030271835.1 |
| :-: | :-: | :-: | :-: | :-: | :-: |
| (1) GCA_000013265.1 | 1.0             | 0.3410          | 0.0024          | 0.3071          | 0.0024 |
| (2) GCA_000005845.2 |                 | 1.0             | 0.0026          | 0.4365          | 0.0026 |
| (3) GCA_000069965.1 |                 |                 | 1.0             | 0.0023          | 0.0311 |
| (4) GCA_000008865.2 |                 |                 |                 | 1.0             | 0.0023 |
| (5) GCA_030271835.1 |                 |                 |                 |                 | 1.0 |


### Rendu TP 2

La matrice de jaccad doit être racalculée en utilisant **une méthode d'échantillonage** de votre choix parmis celles vues en TD.
En plus des 3 bactéries, vous devrez également comparer les kmers issus de ces 3 organismes :

* [Un Humain](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
* [Une souris](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)
* [Un singe](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_029289425.3/)
