import sys
import pysam
import gzip
from os.path import basename, join, dirname

def main(fastq, genome_sam_r1, genome_sam_r2,
         insertseq_sam_r1, insertseq_sam_r2, class_r1, class_r2,
         class_insertseq_sam_stats, class_genome_sam_stats):

    class_genome_sam_dir = dirname(class_genome_sam_stats)
    class_insertseq_sam_dir = dirname(class_insertseq_sam_stats)

    genome_template = pysam.AlignmentFile(genome_sam_r1, "r")
    insertseq_template = pysam.AlignmentFile(insertseq_sam_r1, "r")

    sample = basename(genome_sam_r1).split('.')[0]
    genome = basename(genome_sam_r1).split('.')[-2]
    inseq = basename(insertseq_sam_r1).split('.')[-2]

    genome_sam_r1 = sam_gen(genome_sam_r1)
    genome_sam_r2 = sam_gen(genome_sam_r2)

    insertseq_sam_r1 = sam_gen(insertseq_sam_r1)
    insertseq_sam_r2 = sam_gen(insertseq_sam_r2)

    class_r1 = class_gen(class_r1)
    class_r2 = class_gen(class_r2)

    fastq = fastq_gen(fastq)

    output_sams = dict()

    genome_read1 = next(genome_sam_r1)
    genome_read2 = next(genome_sam_r2)
    inseq_read1 = next(insertseq_sam_r1)
    inseq_read2 = next(insertseq_sam_r2)

    class_tuple = zip(class_r1, class_r2)

    for class1, class2 in class_tuple:

        c1_name, c1_taxon = class1
        c2_name, c2_taxon = class2

        # Ugly code, fix this later
        c1_name = next(fastq)

        gen1 = genome_read1[0]
        gen2 = genome_read2[0]
        is1 = inseq_read1[0]
        is2 = inseq_read2[0]

        if c1_name == gen1 == gen2 == is1 == is2:
            genome_read1 = next(genome_sam_r1)
            genome_read2 = next(genome_sam_r2)
            inseq_read1 = next(insertseq_sam_r1)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen1 == gen2 == is1:
            genome_read1 = next(genome_sam_r1)
            genome_read2 = next(genome_sam_r2)
            inseq_read1 = next(insertseq_sam_r1)
            continue

        elif c1_name == gen1 == gen2 == is2:
            genome_read1 = next(genome_sam_r1)
            genome_read2 = next(genome_sam_r2)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen1 == is1 == is2:
            genome_read1 = next(genome_sam_r1)
            inseq_read1 = next(insertseq_sam_r1)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen2 == is1 == is2:
            genome_read2 = next(genome_sam_r2)
            inseq_read1 = next(insertseq_sam_r1)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen1 == gen2:
            genome_read1 = next(genome_sam_r1)
            genome_read2 = next(genome_sam_r2)
            continue

        elif c1_name == is1 == is2:
            inseq_read1 = next(insertseq_sam_r1)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen1 == is1:
            genome_read1 = next(genome_sam_r1)
            inseq_read1 = next(insertseq_sam_r1)
            continue

        elif c1_name == gen2 == is2:
            genome_read2 = next(genome_sam_r2)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen1 == is2:

            if not c1_taxon in output_sams:
                output_sams[c1_taxon] = [pysam.AlignmentFile(
                    join(class_genome_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                              genome=genome,
                                                                                              inseq=inseq,
                                                                                              taxon=c1_taxon)),
                    "wb", template=genome_template),
                    pysam.AlignmentFile(
                        join(class_insertseq_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                                     genome=genome,
                                                                                                     inseq=inseq,
                                                                                                     taxon=c1_taxon)),
                        "wb", template=insertseq_template)
                ]

            genread_out = genome_read1[1]
            inseq_read_out = inseq_read2[1]
            genread_out.query_name = genread_out.query_name + ':' + c1_taxon
            inseq_read_out.query_name = inseq_read_out.query_name + ':' + c1_taxon
            output_sams[c1_taxon][0].write(genread_out)
            output_sams[c1_taxon][1].write(inseq_read_out)

            genome_read1 = next(genome_sam_r1)
            inseq_read2 = next(insertseq_sam_r2)
            continue

        elif c1_name == gen2 == is1:

            if not c2_taxon in output_sams:
                output_sams[c2_taxon] = [pysam.AlignmentFile(
                    join(class_genome_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                              genome=genome,
                                                                                              inseq=inseq,
                                                                                              taxon=c2_taxon)),
                    "wb", template=genome_template),
                    pysam.AlignmentFile(
                        join(class_insertseq_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                                     genome=genome,
                                                                                                     inseq=inseq,
                                                                                                     taxon=c2_taxon)),
                        "wb", template=insertseq_template)
                ]

            genread_out = genome_read2[1]
            inseq_read_out = inseq_read1[1]
            genread_out.query_name = genread_out.query_name + ':' + c2_taxon
            inseq_read_out.query_name = inseq_read_out.query_name + ':' + c2_taxon
            output_sams[c2_taxon][0].write(genread_out)
            output_sams[c2_taxon][1].write(inseq_read_out)

            genome_read2 = next(genome_sam_r2)
            inseq_read1 = next(insertseq_sam_r1)
            continue

        elif c1_name == gen1:
            genome_read1 = next(genome_sam_r1)

        elif c1_name == gen2:
            genome_read2 = next(genome_sam_r2)

        elif c1_name == is1:
            inseq_read1 = next(insertseq_sam_r1)

        elif c1_name == is2:
            inseq_read2 = next(insertseq_sam_r2)

        else:
            continue

    with open(join(class_genome_sam_dir, "{sample}.{genome}.{inseq}.stats".format(sample=sample,
                                                                                  genome=genome,
                                                                                  inseq=inseq)), 'w') as stats_out1:
        stats_out1.write("STATS GO HERE")

    with open(join(class_insertseq_sam_dir, "{sample}.{genome}.{inseq}.stats".format(sample=sample,
                                                                                     genome=genome,
                                                                                     inseq=inseq)), 'w') as stats_out2:
        stats_out2.write("STATS GO HERE")


def is_mapped(read):
    return not read.is_unmapped


def sam_gen(sam_path):
    for read in pysam.AlignmentFile(sam_path, "r"):
        yield (read.query_name, read)
    while True:
        yield (-1, -1)


def class_gen(class_path):
    with gzip.open(class_path, 'rt') as infile:
        infile.readline()
        for line in infile:
            line=line.strip().split()
            yield ('tmp', line[0])

def fastq_gen(fastq_path):
    with gzip.open(fastq_path, 'rt') as infile:
        num = 0
        for line in infile:
            if num % 4 == 0:
                yield line.strip().split()[0].strip('@')
            num += 1


if __name__ == '__main__':
    fastq = sys.argv[1]
    genome_sam_r1 = sys.argv[2]
    genome_sam_r2 = sys.argv[3]
    insertseq_sam_r1 = sys.argv[4]
    insertseq_sam_r2 = sys.argv[5]
    class_r1 = sys.argv[6]
    class_r2 = sys.argv[7]
    class_genome_sam_stats = sys.argv[8]
    class_insertseq_sam_stats = sys.argv[9]


    main(fastq, genome_sam_r1, genome_sam_r2, insertseq_sam_r1,
         insertseq_sam_r2, class_r1, class_r2,
         class_insertseq_sam_stats, class_genome_sam_stats)