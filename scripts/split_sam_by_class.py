import sys
import pysam
from os.path import basename, join, dirname

def main(genome_sam_r1, genome_sam_r2, insertseq_sam_r1, insertseq_sam_r2, class_r1, class_r2,
         class_insertseq_sam_stats, class_genome_sam_stats):

    class_genome_sam_dir = dirname(class_genome_sam_stats)
    class_insertseq_sam_dir = dirname(class_insertseq_sam_stats)


    sample = basename(genome_sam_r1).split('.')[0]
    genome = basename(genome_sam_r1).split('.')[-2]
    inseq = basename(insertseq_sam_r1).split('.')[-2]

    genome_sam_r1 = pysam.AlignmentFile(genome_sam_r1, "r")
    genome_sam_r2 = pysam.AlignmentFile(genome_sam_r2, "r")

    insertseq_sam_r1 = pysam.AlignmentFile(insertseq_sam_r1, "r")
    insertseq_sam_r2 = pysam.AlignmentFile(insertseq_sam_r2, "r")

    class_r1 = class_gen(class_r1)
    class_r2 = class_gen(class_r2)

    output_sams = dict()

    zip_tuple = zip(genome_sam_r1, genome_sam_r2, insertseq_sam_r1, insertseq_sam_r2, class_r1, class_r2)

    for genome_read1, genome_read2, inseq_read1, inseq_read2, class1, class2 in zip_tuple:

        if not (genome_read1.query_name == genome_read2.query_name == class1[0] == class2[0]):
            print("ERROR: The fastq and class files you provided are not ordered properly. Aborting")
            print("CASE:", sample, genome, inseq)
            sys.exit()

        if genome_read1.is_unmapped and inseq_read2.is_unmapped and genome_read2.is_unmapped and inseq_read1.is_unmapped:
            continue

        elif is_mapped(genome_read1) and is_mapped(inseq_read2):
            if not class1[1] in output_sams:
                output_sams[class1[1]] = [pysam.AlignmentFile(
                    join(class_genome_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                              genome=genome,
                                                                                              inseq=inseq,
                                                                                              taxon=class1[1])),
                    "wb", template=genome_sam_r1),
                    pysam.AlignmentFile(
                        join(class_insertseq_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                                     genome=genome,
                                                                                                     inseq=inseq,
                                                                                                     taxon=class1[1])),
                        "wb", template=insertseq_sam_r2)
                ]

            genome_read1.query_name = genome_read1.query_name + ':' + class1[1]
            output_sams[class1[1]][0].write(genome_read1)
            output_sams[class1[1]][1].write(inseq_read2)

        elif is_mapped(genome_read2) and is_mapped(inseq_read1):
            if not class2[1] in output_sams:
                output_sams[class2[1]] = [pysam.AlignmentFile(
                    join(class_genome_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                              genome=genome,
                                                                                              inseq=inseq,
                                                                                              taxon=class2[1])),
                    "wb", template=genome_sam_r2),
                    pysam.AlignmentFile(
                        join(class_insertseq_sam_dir, "{sample}.{genome}.{inseq}.{taxon}.bam".format(sample=sample,
                                                                                                     genome=genome,
                                                                                                     inseq=inseq,
                                                                                                     taxon=class2[1])),
                        "wb", template=insertseq_sam_r1)
                ]
            genome_read2.query_name = genome_read2.query_name + ':' +  class2[1]
            output_sams[class2[1]][0].write(genome_read2)
            output_sams[class2[1]][1].write(inseq_read1)


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

def class_gen(class_path):
    with open(class_path) as infile:
        infile.readline()
        for line in infile:
            line=line.strip().split()
            yield(line[0].strip('@'), line[2])


if __name__ == '__main__':

    genome_sam_r1 = sys.argv[1]
    genome_sam_r2 = sys.argv[2]
    insertseq_sam_r1 = sys.argv[3]
    insertseq_sam_r2 = sys.argv[4]
    class_r1 = sys.argv[5]
    class_r2 = sys.argv[6]
    class_genome_sam_stats = sys.argv[7]
    class_insertseq_sam_stats = sys.argv[8]


    main(genome_sam_r1, genome_sam_r2, insertseq_sam_r1, insertseq_sam_r2, class_r1, class_r2,
         class_insertseq_sam_stats, class_genome_sam_stats)