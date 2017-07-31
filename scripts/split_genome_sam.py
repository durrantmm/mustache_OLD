import pysam
import sys
from os.path import basename
from operator import itemgetter

def main():
    sam_in = sys.argv[1]
    sams_out = sys.argv[2:]

    samfile = pysam.AlignmentFile(sam_in, "r")

    split_sam(samfile, sams_out)

def split_sam(samfile, sams_out):


    genomes = get_genome_names(sams_out)

    out_sams = get_output_sams(sams_out, samfile)

    dup_reads = duplicate_read_container(next(samfile), genomes)

    for read in samfile:

        name = read.query_name
        contig = read.reference_name
        genome = contig.split('.')[0]

        if name == dup_reads['read_name']:
            dup_reads[genome].append(read)

        else:
            best_reads = get_best_reads(dup_reads, genomes)

            for genome in best_reads.keys():
                if best_reads[genome]:
                    out_sams[genome].write(best_reads[genome])

            dup_reads = duplicate_read_container(read, genomes)


def duplicate_read_container(first_read, genomes):
    out = {genome:[] for genome in genomes}

    name = first_read.query_name
    contig = first_read.reference_name
    genome = contig.split('.')[0]
    out['read_name'] = name

    out[genome].append(first_read)
    return out


def get_output_sams(sams_out, samfile):
    out = dict()

    for sam in sams_out:
        genome = basename(sam).split('.')[-2]
        out[genome] = pysam.AlignmentFile(sam, "wb", template=samfile)

    return out


def get_genome_names(sams_out):
    return [basename(sam).split('.')[-2] for sam in sams_out]


def get_best_reads(dup_reads, genomes):
    out = {genome:None for genome in genomes}

    for genome in genomes:
        scores = [(read.get_tag('AS'), read) for read in dup_reads[genome]]
        scores.sort(key=itemgetter(0), reverse=True)

        if len(scores) > 1:

            if scores[0][0] == scores[1][0]: # Two reads mapping to different loci have same alignment score.
                out[genome] = None
            else:
                out[genome] = scores[0][1]

        elif len(scores) == 1:
            out[genome] = scores[0][1]

        else:
            out[genome] = None


    return out




if __name__ == '__main__':

    main()