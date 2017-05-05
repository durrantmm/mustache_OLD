import sys
import subprocess
from os.path import join, dirname, basename, abspath
from glob import glob


def main(genomes_dir, sample, genome, inseq, taxonomy_dir, outbam):
    taxon_nodes, taxon_ranks = get_taxon_nodes_ranks(join(taxonomy_dir, 'nodes.dmp'), join(taxonomy_dir, 'merged.dmp'))
    inbams = glob(join(genomes_dir, '{sample}.{genome}.{inseq}.*.bam'.format(sample=sample, genome=genome, inseq=inseq)))
    genome_hierarchy = get_taxon_hierarchy_list(genome, taxon_nodes)
    taxa = get_taxa_from_bams(inbams)
    related_taxa = filter_taxa_by_hierarchy(taxa, genome, genome_hierarchy, taxon_nodes)

    bam_files = get_bam_files_to_merge(genomes_dir, sample, genome, inseq, related_taxa)
    samtools_merge(outbam, bam_files)


def samtools_merge(out_bam, bam_files):
    args = "samtools merge {outbam} {inbams}".format(outbam=out_bam, inbams=" ".join(bam_files))
    subprocess.call(args.split())


def get_bam_files_to_merge(genomes_dir, sample, genome, inseq, taxa):
    bam_paths = []
    for taxon in taxa:
        bam_path = join(genomes_dir, '{sample}.{genome}.{inseq}.{taxon}.bam'.format(sample=sample,
                                                                                    genome=genome,
                                                                                    inseq=inseq,
                                                                                    taxon=taxon))
        bam_paths.append(bam_path)
    return bam_paths


def filter_taxa_by_hierarchy(taxa, start_taxon, hierarchy, taxon_nodes):
    out_taxa = []
    for taxon in taxa:
        hierarchy2 = get_taxon_hierarchy_list(taxon, taxon_nodes)
        if is_parent_child(start_taxon, hierarchy, taxon, hierarchy2) or start_taxon == taxon:
            out_taxa.append(taxon)

    return out_taxa


def is_parent_child(start_taxon1, hierarchy1, start_taxon2, hierarchy2):
    if (start_taxon1 not in hierarchy2) and (start_taxon2 in hierarchy1):
        return True
    if (start_taxon2 not in hierarchy1) and (start_taxon1 in hierarchy2):
        return True
    return False


def get_taxon_hierarchy_list(taxon_id, taxon_nodes_dict):
    hierarchy = [taxon_id]

    while taxon_id != '1' and taxon_id != '0':
        taxon_id = taxon_nodes_dict[taxon_id]
        hierarchy.append(taxon_id)

    return hierarchy


def get_taxa_from_bams(inbams):
    taxa = []
    for bamfile in inbams:
        bamfile = basename(bamfile)
        taxa.append(bamfile.split('.')[-2])
    return taxa


def get_taxon_nodes_ranks(nodes_path, merged_path):

    taxon_nodes_dict = {}
    taxon_ranks_dict = {}
    with open(nodes_path) as nodes_in:
        for line in nodes_in:
            line = line.strip().split("|")
            id = line[0].strip()
            rank = line[2].strip()
            parent_id = line[1].strip()
            taxon_nodes_dict[id] = parent_id
            taxon_ranks_dict[id] = rank

    with open(merged_path) as merged_in:
        for line in merged_in:
            line = line.strip().split("|")
            orig_id = line[0].strip()
            merged_id = line[1].strip()
            taxon_nodes_dict[orig_id] = merged_id
            taxon_ranks_dict[orig_id] = taxon_ranks_dict[merged_id]

    return taxon_nodes_dict, taxon_ranks_dict


def get_taxon_hierarchy_list(taxon_id, taxon_nodes_dict):
    hierarchy = [taxon_id]

    while taxon_id != '1' and taxon_id != '0':
        taxon_id = taxon_nodes_dict[taxon_id]
        hierarchy.append(taxon_id)

    return hierarchy


if __name__ == '__main__':

    stats_file = sys.argv[1]
    outbam = sys.argv[2]
    taxonomy_dir= abspath(sys.argv[3])

    genomes_dir = dirname(stats_file)
    sample = basename(stats_file).split('.')[0]
    genome = basename(stats_file).split('.')[1]
    inseq = basename(stats_file).split('.')[2]


    main(genomes_dir, sample, genome, inseq, taxonomy_dir, outbam)