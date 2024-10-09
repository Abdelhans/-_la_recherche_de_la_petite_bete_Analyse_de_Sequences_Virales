
import multiprocessing as mp
import gzip
import re
from Bio import SeqIO
import parasail
import argparse
import time 
import resource 

def generate_suffix_table(sequence):
    seq = str(sequence) + "$"
    tab_suf = []
    for i in range(len(seq)):
        pos = i
        tab_suf.append(pos)
    tab_suf = sorted(tab_suf, key=lambda x: seq[x+1:])
    return tab_suf

def align_read_to_genome(read_seq, genome, suffix_table, k, min_identity):
    kmers = set([read_seq[i:i+k] for i in range(len(read_seq)-k+1)])
    positions = []
    for kmer in kmers:
        start = 0
        end = len(genome) - 1
        while start <= end:
            mid = (start + end) // 2
            if genome[suffix_table[mid]:suffix_table[mid]+k] == kmer:
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
    results = []
    for pos in positions:
        x = read_seq
        z = genome[pos:pos+len(read_seq)]
        result = parasail.sg_dx_scan_sat(x, z, 10, 1, parasail.dnafull)
        identity = float(result.score) / len(read_seq) * 100
        if identity >= min_identity:
            results.append((pos, result.score, identity))
    return results


def align_reads_to_genome(genome_file, reads_file, output_file, k, min_identity, n_processes):
    # Load genome sequences into memory
    genome_sequences = list(SeqIO.parse(genome_file, "fasta"))
    genome = "".join([str(seq.seq) for seq in genome_sequences])
    
    # Generate suffix table for genome
    suffix_table = generate_suffix_table(genome)
    
    # Use a regex to validate sequences and extract kmers
    pattern = re.compile(r'^[ACGT]+$')
    
    # Open output file for writing
    with open(output_file, "w") as out:
        # Write header line
        out.write("read_id\tposition\tscore\tidentity\tsequence\n")
        
        # Open fastq file for reading
        #Ouvre le fichier reads_file en mode lecture avec le module gzip pour lire des fichiers compressés gzip. Crée un objet fichier nommé f pour lire le contenu du fichier.
        with gzip.open(reads_file, "rt") as f:
            #Crée un objet Pool de multiprocessing avec n_processes processus pour traiter les tâches en parallèle.
            pool = mp.Pool(n_processes)
            results = []
            for record in SeqIO.parse(f, "fastq"):
                # Check that sequence is valid
                #Vérifie si la séquence de l'enregistrement record correspond à la regex stockée dans pattern.
                #Si la séquence correspond, le code suivant sera exécuté.
                if pattern.match(str(record.seq)):
                    #Ajoute une tâche au processus pool pour aligner la séquence de l'enregistrement record sur le génome. Stocke le résultat de la tâche dans la liste results.
                    results.append(pool.apply_async(align_read_to_genome, args=(str(record.seq), genome, suffix_table, k, min_identity)))
            #Ferme le processus pool pour indiquer que toutes les tâches ont été soumises et attend que toutes les tâches soient terminées avant de passer à l'étape suivante.        
            pool.close()
            pool.join()
            #Itère sur chaque résultat stocké dans la liste results. 
            for res in results:
                #Pour chaque résultat, itère sur chaque position, score et identité dans la sortie de la tâche en utilisant la méthode get()
                for pos, score, identity in res.get():
                    #Écrit une ligne dans le fichier out avec les valeurs de record.id, pos, score, identity et record.seq
                    out.write(f"{record.id}\t{pos}\t{score}\t{identity:.2f}\t{record.seq}\n")
    
    # Count number of reads aligned
    with open(output_file, "r") as f:
        num_aligned_reads = sum(1 for line in f) - 1  # subtract 1 for header line
        
    elapsed_time = time.process_time()
    mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    
    print(f"{num_aligned_reads} reads on été alignés en {elapsed_time:.2f} seconds, en utilisant  {mem_usage:.2f}")
    
def main():
    
    parser = argparse.ArgumentParser(description='Align reads to genome and write results to output file.')
    parser.add_argument("-g","--genome", type=str, help='path to genome fasta file')
    parser.add_argument("-r","--read", type=str, help='path to input reads fastq.gz file')
    parser.add_argument("-o", "--outfile",type=str, help='path to output results file')
    parser.add_argument("-k", "--kmer_size",type=int, default=11, help='k-mer size (default: 11)')
    parser.add_argument("-t","--threshold", type=float, default=90, help='minimum alignment identity percentage (default: 90)')
    parser.add_argument("-n", "--nbr_pocesseur", type=int, default=1, help='number of processes to use for alignment (default: 1)')

    args = parser.parse_args()

    align_reads_to_genome(args.genome, args.read, args.outfile, args.kmer_size, args.threshold, args.nbr_pocesseur)

if __name__ == '__main__':
    main()
    start_time = time.time()