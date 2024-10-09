# À la recherche de la petite bête

## Description
Le projet "À la recherche de la petite bête" vise à identifier un pathogène mystérieux dans des échantillons pulmonaires de personnes présentant des maladies respiratoires atypiques. Ce projet utilise des techniques de bioinformatique pour aligner des reads sur plusieurs génomes viraux, notamment ceux de la grippe, des rhinovirus, du VIH et des coronavirus, afin de déterminer lequel est le plus similaire au pathogène présent dans l'échantillon.

## Objectifs
1. Aligner des reads sur des génomes de référence.
2. Identifier le virus correspondant au pathogène dans l'échantillon.
3. Filtrer les reads non pertinents pour améliorer la qualité de l'analyse.
4. Analyser la fréquence des k-mers parmi les reads non filtrés pour identifier d'autres pathogènes potentiels.

## Données disponibles
Les données de séquençage et les séquences de référence peuvent être téléchargées à partir des liens suivants :
- Grippe A : [https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/latest_assembly_versions/GCF_000851145.1_ViralMultiSegProj14892/](#)
- VIH-1 : [https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Human_immunodeficiency_virus_1/latest_assembly_versions/GCF_000864765.1_ViralProj15476/GCF_000864765.1_ViralProj15476_cds_from_genomic.fna.gz](#)
- Rhinovirus : [https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Rhinovirus_A/latest_assembly_versions/GCF_000862245.1_ViralProj15330/GCF_000862245.1_ViralProj15330_cds_from_genomic.fna.gz](#)
- Alpha-coronavirus : [https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Alphacoronavirus_1/latest_assembly_versions/GCF_000856025.1_ViralProj15097/GCF_000856025.1_ViralProj15097_cds_from_genomic.fna.gz](#)

## Méthodologie
1. **Extraction de k-mers** : Les séquences sont analysées pour extraire des k-mers, qui sont ensuite utilisés pour l'alignement.
2. **Alignement des reads** : Utilisation d'une stratégie "seed and extend" pour aligner les reads sur les génomes de référence.
3. **Filtrage des reads** : Les reads avec une faible identité ou de nombreuses erreurs sont filtrés.
4. **Analyse des k-mers** : Identification des k-mers les plus fréquents parmi les reads non filtrés pour déceler d'autres pathogènes potentiels.

## Installation
Pour exécuter le projet, vous aurez besoin de Python et des bibliothèques suivantes :
- Biopython
- parasail
- tqdm

Vous pouvez installer les dépendances avec pip :
```bash
pip install biopython parasail tqdm
```
## Utilisation
```bash
python aligner.py --genome <chemin/vers/genome.fasta> --reads <chemin/vers/reads.fastq.gz> --out <chemin/vers/sortie.txt>
```

## Paramètres
--genome : Chemin vers le fichier FASTA contenant le génome de référence.

--reads : Chemin vers le fichier FASTQ.gz contenant les reads à aligner.

--out : Fichier de sortie pour les résultats de l'alignement.

-k : (optionnel) Taille des k-mers utilisés (par défaut : 11).


### Conseils
- Remplacez les liens dans la section "Données disponibles" par les URL réelles où les données peuvent être téléchargées.
- Ajoutez des instructions spécifiques sur la manière de lancer les scripts ou d'utiliser le code si nécessaire.
- Vous pouvez ajouter des sections supplémentaires si vous avez des informations ou des fonctionnalités supplémentaires à inclure. 

N'hésitez pas à personnaliser ce modèle en fonction de vos besoins et du style de votre projet !

