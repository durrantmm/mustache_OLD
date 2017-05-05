import sys
import pysam
from os.path import basename, join
import os
from collections import defaultdict


def main(bam_path, out_path, depth_out_path, lower_perc_depth_cutoff=0.01, lower_depth_cutoff=10, window_size=1000):

    sample = basename(bam_path).split('.')[0]
    genome_taxon = basename(bam_path).split('.')[1]
    inseq = basename(bam_path).split('.')[2]

    bam = pysam.AlignmentFile(bam_path, "rb")

    genome = generate_genome(bam)
    genome, reads = add_reads_to_genome(bam, genome)
    max_depth = get_maximum_depth(genome, reads)

    peaks = get_peaks_scanning_window(genome, reads, max_depth, lower_perc_depth_cutoff, lower_depth_cutoff, window_size)
    peaks = merge_peaks(peaks, genome)

    bam = pysam.AlignmentFile(bam_path, "rb")
    read_stats, read_taxa = get_read_stats(peaks, bam, genome)

    write_depths(peaks, genome, sample, genome_taxon, inseq, depth_out_path)
    write_results(sample, genome_taxon, inseq, read_stats, read_taxa, out_path)


def write_depths(peaks, genome, sample, genome_taxon, inseq, depth_out_path):
    os.makedirs(depth_out_path, exist_ok=True)
    for contig in peaks:
        for start, end in peaks[contig]:
            out_path = join(depth_out_path, "{contig}.{peakstart}.{peakend}.tsv".format(contig=contig,
                                                                                        peakstart=start,
                                                                                        peakend=end))
            with open(out_path, 'w') as out:
                counts = genome[contig]["COUNTS"]
                for_counts = genome[contig]["FOR_COUNTS"]
                rev_counts = genome[contig]["REV_COUNTS"]
                valid_counts = genome[contig]["VALID_COUNTS"]

                zipper = zip(counts[start:end], valid_counts[start:end], for_counts[start:end], rev_counts[start:end])
                i=start

                header = ['CONTIG', 'POS', 'TOTAL_DEPTH', 'VALID_DEPTH', 'FORWARD_DEPTH', 'REVERSE_DEPTH']
                out.write('\t'.join(header)+'\n')
                for tot_count, valid_count, for_count, rev_count in zipper:
                    outline = [contig, i, tot_count, valid_count, for_count, rev_count]
                    outline = map(str, outline)
                    out.write("\t".join(outline)+'\n')
                    i += 1


def write_results(sample, genome_taxon, inseq, read_stats, taxa, out_path):
    outfile = open(out_path, 'w')
    header = ['Sample', 'GenomeTaxon', 'InsertionSeq', 'Contig', 'PeakStart', 'PeakEnd', 'MaxDepth', 'ReadCount', genome_taxon]
    header = header + [taxon for taxon in taxa if taxon != genome_taxon]
    outfile.write("\t".join(header)+'\n')
    for contig in read_stats.keys():
        for peak in read_stats[contig]:
            stats = read_stats[contig][peak]
            outline = [sample, genome_taxon, inseq, contig, peak[0], peak[1], stats['MAX_DEPTH'], stats['READ_COUNT']]

            if genome_taxon in stats['READ_TAXA']:
                outline.append(stats['READ_TAXA'][genome_taxon])
            else:
                outline.append(0)

            for taxon in taxa:
                if taxon != genome_taxon:
                    if taxon in stats['READ_TAXA']:
                        outline.append(stats['READ_TAXA'][taxon])
                    else:
                        outline.append(0)
                else:
                    continue


            outline = map(str, outline)
            outfile.write("\t".join(outline)+'\n')




def get_read_stats(peaks, bam, genome):

    all_read_taxa = set()
    read_stats = defaultdict(lambda: defaultdict(dict))
    for contig in peaks.keys():

        counts = genome[contig]['COUNTS']

        for peak in peaks[contig]:
            start, end = peak
            max_depth = max(counts[start: end])
            read_stats[contig][peak]["MAX_DEPTH"] = max_depth
            read_stats[contig][peak]["READ_COUNT"] = len(list(bam.fetch(contig, start, end)))

            read_taxa = get_read_taxa(list(bam.fetch(contig, start, end)))

            all_read_taxa = all_read_taxa.union(read_taxa.keys())
            read_stats[contig][peak]["READ_TAXA"] = read_taxa

    return dict(read_stats), all_read_taxa


def get_read_taxa(bam_fetch):
    taxa_counts = defaultdict(int)
    for read in bam_fetch:
        taxon = read.query_name.split(':')[-1]
        taxa_counts[taxon] += 1
    return dict(taxa_counts)


def merge_peaks(peaks, genome):
    peaks_out = defaultdict(list)

    for contig in peaks.keys():
        top_peak = None
        for peak in peaks[contig]:

            if not top_peak:
                top_peak = peak
            else:

                if is_overlap(peak, top_peak):
                    counts = genome[contig]['COUNTS']
                    if sum(counts[peak[0]:peak[1]]) > sum(counts[top_peak[0]:top_peak[1]]):
                        top_peak = peak
                    continue

                else:
                    top_peak = trim_peak_by_depth(contig, top_peak, genome)
                    peaks_out[contig].append(top_peak)
                    top_peak = peak

        top_peak = trim_peak_by_depth(contig, top_peak, genome)
        peaks_out[contig].append(top_peak)

    return peaks_out


def trim_peak_by_depth(contig, peak, genome, min_depth = 1):
    counts = genome[contig]['COUNTS']
    peak_start = peak[0]
    peak_end = peak[1]

    while counts[peak_start] < min_depth:
        peak_start += 1

    while counts[peak_end] < min_depth:
        peak_end -= 1

    return (peak_start, peak_end)


def is_overlap(peak1, peak2):
    if peak1[0] <= peak2[1] and peak2[0] <= peak1[1]:
        return True
    else:
        return False


def get_peaks_scanning_window(genome, reads, max_depth, lower_perc_depth_cutoff, lower_depth_cutoff, window_size):
    peaks = defaultdict(list)

    for contig in reads.keys():

        for read in reads[contig]:

            counts = genome[contig]['COUNTS']

            true_start = read.reference_start - read.query_alignment_start

            count_window = counts[true_start: true_start+window_size]

            if max(count_window) > lower_depth_cutoff and max(count_window) > lower_perc_depth_cutoff*max_depth:
                peaks[contig].append((true_start, true_start+window_size))
    return peaks


def get_maximum_depth(genome, reads):

    max_depth = 0
    for contig in genome.keys():
        counts = genome[contig]['COUNTS']
        contig_max_depth = max(counts)
        if contig_max_depth > max_depth:
            max_depth = contig_max_depth

    return max_depth


def add_reads_to_genome(bam, genome):
    reads = defaultdict(list)
    read_count = 0

    contig = ""
    for read in bam:
        new_contig = read.reference_name

        if new_contig != contig:
            read_count = 0

        contig = new_contig

        true_start = read.reference_start - read.query_alignment_start
        true_end = true_start + read.query_length

        for pos in range(true_start, true_end):

            genome[contig]['READS'][pos].add(read_count)
            genome[contig]['COUNTS'][pos] += 1
            if not read.is_reverse:
                if pos >= read.reference_start and pos < read.reference_end:
                    genome[contig]["FOR_COUNTS"][pos] += 1
                    genome[contig]["VALID_COUNTS"][pos] += 1

            else:
                if pos >= read.reference_start and pos < read.reference_end:
                    genome[contig]["REV_COUNTS"][pos] += 1
                    genome[contig]["VALID_COUNTS"][pos] += 1

        reads[contig].append(read)
        read_count += 1

    return genome, reads



def generate_genome(bam):
    genome = defaultdict(dict)
    for contig in bam.header['SQ']:
        genome[contig['SN']]["READS"] = [set()] * contig['LN']
        genome[contig['SN']]["COUNTS"] = [0] * contig['LN']
        genome[contig['SN']]["VALID_COUNTS"] = [0] * contig['LN']
        genome[contig['SN']]["FOR_COUNTS"] = [0] * contig['LN']
        genome[contig['SN']]["REV_COUNTS"] = [0] * contig['LN']
    return genome


if __name__ == '__main__':
    bam_path = sys.argv[1]
    out_path = sys.argv[2]
    depth_out_path = sys.argv[3]

    main(bam_path, out_path, depth_out_path)