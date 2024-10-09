import gzip
import re
from Bio import SeqIO
import parasail
import argparse
from tqdm import tqdm
import resource
import time
import sys
import doctest


def extract_kmers_from_fastq_file(filename, k):
    """
    Cette fonction prend en entrée un fichier fastq.gz et un entier k, et retourne un set contenant tous les k-mers
    uniques extraits des séquences du fichier. Si une séquence contient des caractères autres que A, C, G ou T, elle est
    ignorée. Les k-mers sont extraits en parcourant toutes les sous-séquences de longueur k dans chaque séquence.
    
    Args:
        filename (str): le nom du fichier fastq.gz à lire
        k (int): la longueur des k-mers à extraire
        
    Returns:
        set: un set contenant tous les k-mers uniques extraits des séquences du fichier
    
    Examples:
        # On teste la fonction sur un petit fichier fastq.gz contenant deux séquences
        >>> test_file = 'test.fastq.gz'
        >>> with gzip.open(test_file, 'wt') as f:
        ...     f.write('@seq1\nACGTACGTACGT\n+\n????????????\n@seq2\nACGTNGTACGT\n+\n????????????\n')
        >>> extract_kmers_from_fastq_file(test_file, 3)
        {'CGT', 'ACG', 'TAC', 'GTN', 'CGA', 'GTA', 'ACN', 'CTA', 'NGT'}
        
        # On teste la fonction sur un fichier fastq.gz qui n'existe pas
        >>> extract_kmers_from_fastq_file('file_not_found.fastq.gz', 3)
        Traceback (most recent call last):
            ...
        FileNotFoundError: [Errno 2] No such file or directory: 'file_not_found.fastq.gz'
        
        # On teste la fonction avec une valeur de k invalide
        >>> extract_kmers_from_fastq_file(test_file, 0)
        Traceback (most recent call last):
            ...
        ValueError: k doit être un entier positif non nul
    """
    if not isinstance(k, int) or k <= 0:
        raise ValueError('k doit être un entier positif non nul')
    # On ouvre le fichier fastq.gz en mode lecture avec gzip
    with gzip.open(filename, "rt") as f:
        # On initialise un set pour stocker les k-mers uniques
        kmers = set()
        # On utilise une expression régulière pour valider les séquences et extraire les k-mers
        pattern = re.compile(r'^[ACGT]+$')
        # On parcourt le fichier ligne par ligne
        for line_num, line in enumerate(f):
            # On traite les lignes correspondant aux séquences
            if line_num % 4 == 1:
                # On retire les caractères de saut de ligne et on transforme la séquence en majuscules
                sequence = line.strip().upper()
                # On vérifie que la séquence ne contient que des caractères valides (A, C, G, T)
                if pattern.match(sequence):
                    # On ajoute tous les k-mers de la séquence au set
                    kmers.update([sequence[i:i+k] for i in range(len(sequence)-k+1)])
    # On retourne la liste des k-mers
    return kmers
def extract_sequences_from_fasta_file(filename):
    """
    Extracts sequences from a FASTA file and returns them as a list of strings.

    Args:
    filename (str): the name of the FASTA file to extract sequences from.

    Returns:
    list: a list of strings, where each string is a sequence from the FASTA file.

    Raises:
    ValueError: if the provided filename does not have a .fasta or .fa extension.

    Examples:
    >>> extract_sequences_from_fasta_file("example.fasta")
    ['ATGCTAGCTAGCTAG', 'CGATCGATCGAT', 'ATCGATCGATCGATCGAT']
    """

    # Check that the file has a valid extension
    if not filename.endswith(".fasta") and not filename.endswith(".fa"):
        raise ValueError("Invalid file format. Expected .fasta or .fa extension.")
    sequences = []
    with open(filename, "r") as handle:
        # Ouvre le fichier `filename` en mode lecture (r) et crée un objet `handle` pour y accéder.
        # Utilise la structure de bloc `with` pour s'assurer que le fichier est fermé automatiquement à la fin.
        for record in SeqIO.parse(handle, "fasta"):
            # Parcourt les enregistrements dans le fichier FASTA en utilisant la fonction `SeqIO.parse()` de BioPython.
            # Chaque enregistrement est stocké dans l'objet `record`.
            sequences.append(str(record.seq))
            # Extrait la séquence de l'enregistrement et la stocke sous forme de chaîne de caractères dans la liste `sequences`.
    return sequences

def generate_suffix_table(sequence):
    """
    Cette fonction prend en entrée une séquence et retourne une table des suffixes trié dans l'ordre lexicologique.

    Args:
        sequence (str): La séquence d'ADN ou d'ARN à partir de laquelle la table des suffixes est généré.

    Returns:
        list: Un tableau d'entiers représentant les positions des suffixes trié lexicologiquement.

    Examples:
        >>> generate_suffix_table("ATCGGACT")
        [8, 0, 5, 2, 1, 6, 3, 7, 4]
        >>> generate_suffix_table("AAACGTAAGCTAGCTTACG$")
        [19, 17, 14, 9, 15, 6, 11, 2, 8, 1, 13, 18, 7, 16, 5, 10, 12, 4, 3, 0]
    """
    seq = str(sequence) + "$" # Ajoute un symbole $ à la fin de la séquence pour marquer la fin.
    tab_suf = []
    for i in range(len(seq)):
        pos = i
        tab_suf.append(pos) # Ajoute la position de chaque suffixe dans la table.
    tab_suf = sorted(tab_suf, key=lambda x: seq[x+1:]) # Trie la table des suffixes en utilisant la fonction lambda selon l'ordre suffixe lexicologique.
    return tab_suf


def search_kmers_in_genome(fastq_file, fasta_file, k):
    """
    Search for all occurrences of kmers of length k in a genome sequence and store their positions.
    
    Args:
    - fastq_file (str): The path to the FASTQ file containing the genome sequence.
    - fasta_file (str): The path to the FASTA file containing the kmers to search for.
    - k (int): The length of the kmers to search for.
    
    Returns:
    - positions (dict): A dictionary containing the positions of each kmer in the genome sequence. 
    The keys of the dictionary are the kmers, and the values are lists of integers representing 
    the positions of the kmer in the genome sequence.
    
    Example:
    >>> genome = "ACGTACGTACGT"
    >>> kmers = ["ACGT", "CGTA", "GTAC"]
    >>> suffix_table = [0, 4, 8]
    >>> search_kmers_in_genome(genome, kmers, suffix_table, 4)
    {'ACGT': [0, 4, 8], 'CGTA': [1, 5], 'GTAC': [2, 6]}
    """
    # Extract kmers from fastq file
    kmers = extract_kmers_from_fastq_file(fastq_file, k)
    
    # Extract sequences from fasta file
    sequences = extract_sequences_from_fasta_file(fasta_file)
    genome = ''.join(sequences)
    
    # Generate suffix table for genome
    suffix_table = generate_suffix_table(genome)
    
    # Search kmers in genome and store positions
    # Create an empty dictionary to store the positions of each kmer
    positions = {}
    # For each kmer in the list of kmers to search for
    for kmer in kmers:
        # Set the start and end positions of the genome sequence to search within
        start = 0
        end = len(genome) - 1
        # Use binary search to find the position of the first occurrence of the kmer in the genome sequence
        while start <= end:
            mid = (start + end) // 2
            if genome[suffix_table[mid]:suffix_table[mid]+k] == kmer:
                # If the kmer is found, expand to left and right to find all occurrences
                left = mid - 1
                while left >= 0 and genome[suffix_table[left]:suffix_table[left]+k] == kmer:
                    left -= 1
                right = mid + 1
                while right < len(genome) and genome[suffix_table[right]:suffix_table[right]+k] == kmer:
                    right += 1
                # Store the positions of the kmer in the dictionary
                positions[kmer] = [suffix_table[i] for i in range(left+1, right)]
                break
            elif genome[suffix_table[mid]:suffix_table[mid]+k] < kmer:
                start = mid + 1
            else:
                end = mid - 1
    return positions


def align_reads_to_genome(genome_file, reads_file, output_file, k, threshold):
    # Load genome sequences into memory
    genome_sequences = list(SeqIO.parse(genome_file, "fasta"))
    genome = "".join([str(seq.seq) for seq in genome_sequences])
    
    # Generate suffix table for genome
    suffix_table = generate_suffix_table(genome)
    
    # Open output file for writing
    with open(output_file, "w") as out:
        # Write header line
        out.write("read_id\t-----position\tscore\tidentity\tsequence\n")
        
        # Open fastq file for reading
        with gzip.open(reads_file, "rt") as f:
            # Use a regex to validate sequences and extract kmers
            pattern = re.compile(r'^[ACGT]+$')
            
            # Iterate over fastq records with tqdm progress bar
            for record in tqdm(SeqIO.parse(f, "fastq"), total=10000):
                # Check that sequence is valid
                if pattern.match(str(record.seq)):
                    # Extract kmers from read sequence
                    kmers = set([str(record.seq)[i:i+k] for i in range(len(record.seq)-k+1)])
                    
                    # Search kmers in genome and store positions
                    positions = []
                    for kmer in kmers:
                        start = 0
                        end = len(genome) - 1
                        while start <= end:
                            mid = (start + end) // 2
                            if genome[suffix_table[mid]:suffix_table[mid]+k] == kmer:
                                # Found kmer, expand to left and right to find all occurrences
                                left = mid - 1
                                while left >= 0 and genome[suffix_table[left]:suffix_table[left]+k] == kmer:
                                    left -= 1
                                right = mid + 1
                                while right < len(genome) and genome[suffix_table[right]:suffix_table[right]+k] == kmer:
                                    right += 1
                                positions += [suffix_table[i] for i in range(left+1, right)]
                                break
                            elif genome[suffix_table[mid]:suffix_table[mid]+k] < kmer:
                                start = mid + 1
                            else:
                                end = mid - 1
                    
                    # Align read to genome for each position found
                    # Pour chaque position où un kmer a été trouvé dans la séquence génomique
                    for pos in positions:
                        # Convertir la séquence de l'objet record en chaîne de caractères
                        x = str(record.seq)
                        # Extraire une sous-chaîne de la séquence génomique à partir de la position actuelle
                        # jusqu'à la longueur de la séquence de l'objet record
                        z = genome[pos:pos+len(record.seq)]
                        # Utilise l'algorithme d'alignement de séquences Parasail pour aligner les séquences
                        # de l'objet record et de la sous-chaîne de la séquence génomique
                        result = parasail.sg_dx_scan_sat(x, z, 10, 1, parasail.dnafull)
                        # Calcule du pourcentage d'identité entre les deux séquences alignées
                        identity = float(result.score) / len(record.seq) * 100
                        # Si le pourcentage d'identité est supérieur ou égal au seuil spécifié
                        if identity >= threshold:
                            # On ecrit les informations de l'alignement dans le fichier de sortie au format txt
                            out.write(f"{record.id}\t{pos}\t{result.score}\t{identity:.2f}\t{record.seq}\n")
                   

    # Count number of reads aligned
    with open(output_file, "r") as f:
        # Compter le nombre total d'alignements effectués en comptant le nombre de lignes dans le fichier
        # de sortie, en soustrayant 1 pour l'en-tête
        num_aligned_reads = sum(1 for line in f) - 1  # subtract 1 for header line
    # Mesurer le temps écoulé depuis le début de l'exécution du script en utilisant la fonction process_time
    elapsed_time = time.process_time()
    # Mesurer l'utilisation maximale de la mémoire vive par le processus actuel, en utilisant la fonction
    # getrusage de la bibliothèque resource, et diviser par 1024 pour obtenir le résultat en kilo-octets
    mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    
    print(f"{num_aligned_reads} reads on été alignés en {elapsed_time:.2f} seconds, en utilisant {mem_usage:.2f} ko ")
    
def main():
    parser = argparse.ArgumentParser(description="Align reads to genome")
    parser.add_argument("-g", "--genome", help="Path to genome file in fasta format")
    parser.add_argument("-r", "--read", help="Path to reads file in fastq format")
    parser.add_argument("-o", "--outfile", help="Path to output file")
    parser.add_argument("-k", "--kmer_length", type=int, default=11,
                        help="Length of k-mers to extract from reads (default: 11)")
    parser.add_argument("-t", "--threshold", type=float, default=95,
                        help="Minimum alignment identity threshold (default: 90)")
    args = parser.parse_args()
    align_reads_to_genome(args.genome, args.read, args.outfile, args.kmer_length, args.threshold)

if __name__ == '__main__':
    #doctest.testmod(verbose=True)
    main()
    start_time = time.time()
