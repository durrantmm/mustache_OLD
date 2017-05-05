import sys, os
from Bio import SeqIO
from random import randint, choice

def main(genome_path, insertion_path, number, out_fasta):

    fasta = SeqIO.parse(genome_path, "fasta")

    chrom_lengths = []
    for chrom in fasta:
        chrom_lengths.append(len(chrom.seq))

    chrom_props = [length/float(sum(chrom_lengths)) for length in chrom_lengths]
    chrom_insertion_number = [int(round(prop*number, 0)) for prop in chrom_props]

    if sum(chrom_insertion_number) != number:
        chrom_insertion_number[0] += sum(chrom_insertion_number) - number

    inseq = list(SeqIO.parse(insertion_path, "fasta"))[0].seq
    genome = SeqIO.parse(genome_path, "fasta")

    out_genome = []
    for i, chrom in enumerate(genome):
        seq = chrom.seq
        insertion_sites = list([randint(0, chrom_lengths[i]-1) for j in range(chrom_insertion_number[i])])
        insertion_sites.sort()

        offset = 0
        for site in sorted(insertion_sites):
            print("%s:%d" % (chrom.name, site))
            seq = seq[:site+offset] + inseq + seq[site+offset:]

        chrom.seq = seq
        out_genome.append(chrom)

    with open(out_fasta, 'w') as out:
        SeqIO.write(out_genome, out, "fasta")





if __name__ == '__main__':
    genome_path = sys.argv[1]
    insertion_seq = sys.argv[2]
    number = int(sys.argv[3])
    output_fasta = sys.argv[4]

    main(genome_path, insertion_seq, number, output_fasta)