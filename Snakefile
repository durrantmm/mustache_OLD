configfile: "config.yaml"

import os
from os.path import basename

OUTPUT_DIR = config['output_dir']
FASTQ_DIR = os.path.join(OUTPUT_DIR, config['fastq_dir'])
GENOME_DIR = os.path.join(OUTPUT_DIR, config['genome_dir'])
CLASS_DIR = os.path.join(OUTPUT_DIR, config['class_dir'])
INSERTSEQ_DIR = os.path.join(OUTPUT_DIR, config['insertseq_dir'])

GENOME_SAM_DIR = os.path.join(OUTPUT_DIR, config['genome_sam_dir'])
if not os.path.isdir(GENOME_SAM_DIR):
    os.makedirs(GENOME_SAM_DIR)

INSERTSEQ_SAM_DIR = os.path.join(OUTPUT_DIR, config['insertseq_sam_dir'])
if not os.path.isdir(INSERTSEQ_SAM_DIR):
    os.makedirs(INSERTSEQ_SAM_DIR)

CLASS_INSERTSEQ_SAM_DIR = os.path.join(OUTPUT_DIR, config['class_insertseq_sam_dir'])
if not os.path.isdir(CLASS_INSERTSEQ_SAM_DIR):
    os.makedirs(CLASS_INSERTSEQ_SAM_DIR)

CLASS_GENOME_SAM_DIR = os.path.join(OUTPUT_DIR, config['class_genome_sam_dir'])
if not os.path.isdir(CLASS_GENOME_SAM_DIR):
    os.makedirs(CLASS_GENOME_SAM_DIR)

SAM_MERGED_TAXONOMY = os.path.join(OUTPUT_DIR, config['sam_merged_taxonomy_dir'])
if not os.path.isdir(SAM_MERGED_TAXONOMY):
    os.makedirs(SAM_MERGED_TAXONOMY)

SAM_SORTED_DIR = os.path.join(OUTPUT_DIR, config['sam_sorted_dir'])
if not os.path.isdir(SAM_SORTED_DIR):
    os.makedirs(SAM_SORTED_DIR)

PEAKS_DIR = os.path.join(OUTPUT_DIR, config['peaks_dir'])
if not os.path.isdir(PEAKS_DIR):
    os.makedirs(PEAKS_DIR)

PEAKS_READ_DEPTHS_DIR = os.path.join(OUTPUT_DIR, config['peaks_read_depth_dir'])
if not os.path.isdir(PEAKS_READ_DEPTHS_DIR):
    os.makedirs(PEAKS_READ_DEPTHS_DIR)

PEAKS_PLOTS_DIR = os.path.join(OUTPUT_DIR, config['peak_plots_dir'])
if not os.path.isdir(PEAKS_PLOTS_DIR):
    os.makedirs(PEAKS_PLOTS_DIR)

WC_fastqs = glob_wildcards(os.path.join(FASTQ_DIR, "{sample}.{pair}.fastq.gz"))
WC_genomes = glob_wildcards(os.path.join(GENOME_DIR, "{genome}.fasta"))
WC_inseqs = glob_wildcards(os.path.join(INSERTSEQ_DIR, "{insertseq}.fasta"))


SAMPLES = set(WC_fastqs.sample)
PAIRS = ['R1', 'R2']

GENOMES = WC_genomes.genome
INSERTSEQS = WC_inseqs.insertseq


rule all:
    input:
        expand("{peak_plots_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peak_plots_dir=PEAKS_PLOTS_DIR),
               sample=SAMPLES, insertseq=INSERTSEQS, genome=GENOMES)
    run:
        print("COMPLETED SUCCESSFULLY!")


rule bowtie_build_genome:
    input:
        "{genome_dir}/{{genome}}.fasta".format(genome_dir=GENOME_DIR)
    output:
        "{genome_dir}/{{genome}}.fasta.1.bt2".format(genome_dir=GENOME_DIR),
        "{genome_dir}/{{genome}}.fasta.2.bt2".format(genome_dir=GENOME_DIR),
        "{genome_dir}/{{genome}}.fasta.3.bt2".format(genome_dir=GENOME_DIR),
        "{genome_dir}/{{genome}}.fasta.4.bt2".format(genome_dir=GENOME_DIR),
        "{genome_dir}/{{genome}}.fasta.rev.1.bt2".format(genome_dir=GENOME_DIR),
        "{genome_dir}/{{genome}}.fasta.rev.2.bt2".format(genome_dir=GENOME_DIR)
    shell:
        "bowtie2-build {input} {input}"


rule bowtie_build_insertseq:
    input:
        "{insertseq_dir}/{{insertseq}}.fasta".format(insertseq_dir=INSERTSEQ_DIR)
    output:
        "{insertseq_dir}/{{insertseq}}.fasta.1.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        "{insertseq_dir}/{{insertseq}}.fasta.2.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        "{insertseq_dir}/{{insertseq}}.fasta.3.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        "{insertseq_dir}/{{insertseq}}.fasta.4.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        "{insertseq_dir}/{{insertseq}}.fasta.rev.1.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        "{insertseq_dir}/{{insertseq}}.fasta.rev.2.bt2".format(insertseq_dir=INSERTSEQ_DIR)
    shell:
        "bowtie2-build {input} {input}"


rule bowtie_align_genome:
    input:
        fastq = "{fastq_dir}/{{sample}}.{{pair}}.fastq.gz".format(fastq_dir=FASTQ_DIR),
        genome = "{genome_dir}/{{genome}}.fasta".format(genome_dir=GENOME_DIR),
        i1 = "{genome_dir}/{{genome}}.fasta.1.bt2".format(genome_dir=GENOME_DIR),
        i2 = "{genome_dir}/{{genome}}.fasta.2.bt2".format(genome_dir=GENOME_DIR),
        i3 = "{genome_dir}/{{genome}}.fasta.3.bt2".format(genome_dir=GENOME_DIR),
        i4 = "{genome_dir}/{{genome}}.fasta.4.bt2".format(genome_dir=GENOME_DIR),
        i5 = "{genome_dir}/{{genome}}.fasta.rev.1.bt2".format(genome_dir=GENOME_DIR),
        i6 = "{genome_dir}/{{genome}}.fasta.rev.2.bt2".format(genome_dir=GENOME_DIR)
    output:
        sam = "{genome_sam_dir}/{{sample}}.{{pair}}.{{genome}}.sam".format(genome_sam_dir=GENOME_SAM_DIR)
    shell:
        "bowtie2 -x {input.genome} -U {input.fastq} -S {output.sam} --local --quiet --reorder"


rule bowtie_align_insertseq:
    input:
        fastq = "{fastq_dir}/{{sample}}.{{pair}}.fastq.gz".format(fastq_dir=FASTQ_DIR),
        insertseq = "{insertseq_dir}/{{insertseq}}.fasta".format(insertseq_dir=INSERTSEQ_DIR),
        i1 = "{insertseq_dir}/{{insertseq}}.fasta.1.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        i2 = "{insertseq_dir}/{{insertseq}}.fasta.2.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        i3 = "{insertseq_dir}/{{insertseq}}.fasta.3.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        i4 = "{insertseq_dir}/{{insertseq}}.fasta.4.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        i5 = "{insertseq_dir}/{{insertseq}}.fasta.rev.1.bt2".format(insertseq_dir=INSERTSEQ_DIR),
        i6 = "{insertseq_dir}/{{insertseq}}.fasta.rev.2.bt2".format(insertseq_dir=INSERTSEQ_DIR)
    output:
        sam = "{insertseq_sam_dir}/{{sample}}.{{pair}}.{{insertseq}}.sam".format(insertseq_sam_dir=INSERTSEQ_SAM_DIR)
    shell:
        "bowtie2 -x {input.insertseq} -U {input.fastq} -S {output.sam} --local --quiet --reorder"


rule split_sam_by_class:
    input:
        genome_sam_R1 = "{genome_sam_dir}/{{sample}}.{pair}.{{genome}}.sam".format(genome_sam_dir=GENOME_SAM_DIR, pair=PAIRS[0]),
        genome_sam_R2 = "{genome_sam_dir}/{{sample}}.{pair}.{{genome}}.sam".format(genome_sam_dir=GENOME_SAM_DIR, pair=PAIRS[1]),
        insertseq_sam_R1 = "{insertseq_sam_dir}/{{sample}}.{pair}.{{insertseq}}.sam".format(insertseq_sam_dir=INSERTSEQ_SAM_DIR, pair=PAIRS[0]),
        insertseq_sam_R2 = "{insertseq_sam_dir}/{{sample}}.{pair}.{{insertseq}}.sam".format(insertseq_sam_dir=INSERTSEQ_SAM_DIR, pair=PAIRS[1]),
        class_R1 = "{class_dir}/{{sample}}.{pair}.class.gz".format(class_dir=CLASS_DIR, pair=PAIRS[0]),
        class_R2 = "{class_dir}/{{sample}}.{pair}.class.gz".format(class_dir=CLASS_DIR, pair=PAIRS[1])

    output:
        genome_sam_stats = "{class_genome_sam_dir}/{{sample}}.{{genome}}.{{insertseq}}.stats".format(class_genome_sam_dir=CLASS_GENOME_SAM_DIR),
        insert_sam_stats = "{class_insertseq_sam_dir}/{{sample}}.{{genome}}.{{insertseq}}.stats".format(class_insertseq_sam_dir=CLASS_INSERTSEQ_SAM_DIR)

    shell:
        "python scripts/split_sam_by_class.py "
        "{input.genome_sam_R1} {input.genome_sam_R2} "
        "{input.insertseq_sam_R1} {input.insertseq_sam_R2} "
        "{input.class_R1} {input.class_R2} "
        "{output.genome_sam_stats} {output.insert_sam_stats}"


rule sam_merge_taxonomy:
    input:
        genome_sam_stats = "{class_genome_sam_dir}/{{sample}}.{{genome}}.{{insertseq}}.stats".format(class_genome_sam_dir=CLASS_GENOME_SAM_DIR),
    output:
        "{sam_merged_taxonomy_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam".format(sam_merged_taxonomy_dir=SAM_MERGED_TAXONOMY)
    shell:
        "python scripts/sam_merge_taxonomy.py {input} {output} " + config['taxonomy_db_dir']


rule samtools_sort:
    input:
        "{sam_merged_taxonomy_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam".format(sam_merged_taxonomy_dir=SAM_MERGED_TAXONOMY)
    output:
        "{sam_sorted_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam".format(sam_sorted_dir=SAM_SORTED_DIR)
    shell:
        "samtools sort {input} > {output}"


rule samtools_index:
    input:
        "{sam_sorted_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam".format(sam_sorted_dir=SAM_SORTED_DIR)
    output:
        "{sam_sorted_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam.bai".format(sam_sorted_dir=SAM_SORTED_DIR)
    shell:
        "samtools index {input}"


rule call_peaks:
    input:
        bam = "{sam_sorted_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam".format(sam_sorted_dir=SAM_SORTED_DIR),
        bai = "{sam_sorted_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam.bai".format(sam_sorted_dir=SAM_SORTED_DIR)
    output:
        "{peaks_dir}/{{sample}}.{{genome}}.{{insertseq}}.peaks.tsv".format(peaks_dir=PEAKS_DIR),
        "{peaks_read_depth_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peaks_read_depth_dir=PEAKS_READ_DEPTHS_DIR)
    shell:
        "python scripts/call_peaks.py {input.bam} {output}"


rule plot_peaks:
    input:
        "{peaks_read_depth_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peaks_read_depth_dir=PEAKS_READ_DEPTHS_DIR)
    output:
        "{peak_plots_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peak_plots_dir=PEAKS_PLOTS_DIR)

    shell:
        "Rscript scripts/create_peak_graphs.R {input} {output}"

