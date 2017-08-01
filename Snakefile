configfile: "config.yaml"

import os
from os.path import basename

OUTPUT_DIR = config['output_dir']
FASTQ_DIR = os.path.join(OUTPUT_DIR, config['fastq_dir'])
GENOME_DIR = os.path.join(OUTPUT_DIR, config['genome_dir'])
CLASS_DIR = os.path.join(OUTPUT_DIR, config['class_dir'])
INSERTSEQ_DIR = os.path.join(OUTPUT_DIR, config['insertseq_dir'])

COMBINED_GENOME_DIR = os.path.join(OUTPUT_DIR, config['combined_genome_dir'])
COMBINED_INSERTSEQ_DIR = os.path.join(OUTPUT_DIR, config['combined_insertseq_dir'])

GENOME_SAM_DIR = os.path.join(OUTPUT_DIR, config['genome_sam_dir'])
INSERTSEQ_SAM_DIR = os.path.join(OUTPUT_DIR, config['insertseq_sam_dir'])

SPLIT_GENOME_SAM_DIR = os.path.join(OUTPUT_DIR, config['split_genome_sam_dir'])
SPLIT_INSERTSEQ_SAM_DIR = os.path.join(OUTPUT_DIR, config['split_insertseq_sam_dir'])

CLASS_INSERTSEQ_SAM_DIR = os.path.join(OUTPUT_DIR, config['class_insertseq_sam_dir'])
CLASS_GENOME_SAM_DIR = os.path.join(OUTPUT_DIR, config['class_genome_sam_dir'])

SAM_MERGED_TAXONOMY = os.path.join(OUTPUT_DIR, config['sam_merged_taxonomy_dir'])
SAM_SORTED_DIR = os.path.join(OUTPUT_DIR, config['sam_sorted_dir'])

PEAKS_DIR = os.path.join(OUTPUT_DIR, config['peaks_dir'])
PEAKS_READ_DEPTHS_DIR = os.path.join(OUTPUT_DIR, config['peaks_read_depth_dir'])
PEAKS_PLOTS_DIR = os.path.join(OUTPUT_DIR, config['peak_plots_dir'])

WC_fastqs = glob_wildcards(os.path.join(FASTQ_DIR, "{sample}.{pair}.fastq.gz"))
WC_genomes = glob_wildcards(os.path.join(GENOME_DIR, "{genome}.fasta"))
WC_inseqs = glob_wildcards(os.path.join(INSERTSEQ_DIR, "{insertseq}.fasta"))


SAMPLES = set(WC_fastqs.sample)
PAIRS = ['R1', 'R2']

GENOMES = WC_genomes.genome
INSERTSEQS = WC_inseqs.insertseq

#print( expand("3.split_genome_sam/{sample}.{pair}.{genome}.sam", sample=SAMPLES, pair=PAIRS, genome=GENOMES))

rule all:
    input:
        #expand("%s/{sample}.{pair}.{genome}.bam" % SPLIT_GENOME_SAM_DIR, sample=SAMPLES, pair=PAIRS, genome=GENOMES),
        #expand("%s/{sample}.{pair}.{insertseq}.bam" % SPLIT_INSERTSEQ_SAM_DIR, sample=SAMPLES, pair=PAIRS, insertseq=INSERTSEQS)
        #expand("{genome_sam_dir}/{{sample}}.{{pair}}.genomes.sam".format(genome_sam_dir=GENOME_SAM_DIR), sample=SAMPLES, pair=PAIRS),
        #expand("{insertseq_sam_dir}/{{sample}}.{{pair}}.insertseqs.sam".format(insertseq_sam_dir=INSERTSEQ_SAM_DIR), sample=SAMPLES, pair=PAIRS)
        expand("{peak_plots_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peak_plots_dir=PEAKS_PLOTS_DIR),
               sample=SAMPLES, insertseq=INSERTSEQS, genome=GENOMES)
    run:
        print("MUSTACHE FINISHED WITH NO EXCEPTIONS!")


rule download_taxonomy_db:
    output:
        config['taxonomy_db_dir']
    shell:
        'python scripts/download_taxonomy_db.py'


rule combine_genomes:
    input:
        expand("{genome_dir}/{{genome}}.fasta".format(genome_dir=GENOME_DIR), genome=GENOMES)
    output:
        "{combined_genome_dir}/genomes.combined.fasta".format(combined_genome_dir=COMBINED_GENOME_DIR)
    shell:
        "python scripts/genome_combine.py {output} {input}"


rule combine_insertseqs:
    input:
        expand("{insertseq_dir}/{{insertseq}}.fasta".format(insertseq_dir=INSERTSEQ_DIR), insertseq=INSERTSEQS)
    output:
        "{combined_insertseq_dir}/insertseqs.combined.fasta".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR)
    shell:
        "python scripts/insertseq_combine.py {output} {input}"


rule bowtie_build_genomes:
    input:
        "{combined_genome_dir}/genomes.combined.fasta".format(combined_genome_dir=COMBINED_GENOME_DIR)
    output:
        "{combined_genome_dir}/genomes.combined.fasta.1.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        "{combined_genome_dir}/genomes.combined.fasta.2.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        "{combined_genome_dir}/genomes.combined.fasta.3.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        "{combined_genome_dir}/genomes.combined.fasta.4.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        "{combined_genome_dir}/genomes.combined.fasta.rev.1.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        "{combined_genome_dir}/genomes.combined.fasta.rev.2.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR)
    shell:
        "bowtie2-build {input} {input}"


rule bowtie_build_insertseqs:
    input:
        "{combined_insertseq_dir}/insertseqs.combined.fasta".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR)
    output:
        "{combined_insertseq_dir}/insertseqs.combined.fasta.1.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        "{combined_insertseq_dir}/insertseqs.combined.fasta.2.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        "{combined_insertseq_dir}/insertseqs.combined.fasta.3.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        "{combined_insertseq_dir}/insertseqs.combined.fasta.4.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        "{combined_insertseq_dir}/insertseqs.combined.fasta.rev.1.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        "{combined_insertseq_dir}/insertseqs.combined.fasta.rev.2.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR)
    shell:
        "bowtie2-build {input} {input}"


rule bowtie_align_to_genomes:
    input:
        fastq = "{fastq_dir}/{{sample}}.{{pair}}.fastq.gz".format(fastq_dir=FASTQ_DIR),
        genomes = "{combined_genome_dir}/genomes.combined.fasta".format(combined_genome_dir=COMBINED_GENOME_DIR),
        in1 = "{combined_genome_dir}/genomes.combined.fasta.1.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        in2 = "{combined_genome_dir}/genomes.combined.fasta.2.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        in3 = "{combined_genome_dir}/genomes.combined.fasta.3.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        in4 = "{combined_genome_dir}/genomes.combined.fasta.4.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        in5 = "{combined_genome_dir}/genomes.combined.fasta.rev.1.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR),
        in6 = "{combined_genome_dir}/genomes.combined.fasta.rev.2.bt2".format(combined_genome_dir=COMBINED_GENOME_DIR)
    output:
        sam = "{genome_sam_dir}/{{sample}}.{{pair}}.genomes.sam".format(genome_sam_dir=GENOME_SAM_DIR)
    threads:
        config['bowtie2_threads']
    shell:
        "bowtie2 -x {input.genomes} -U {input.fastq} -S {output.sam} -p {threads} --local --quiet --reorder --no-unal --all"


rule bowtie_align_to_insertseqs:
    input:
        fastq = "{fastq_dir}/{{sample}}.{{pair}}.fastq.gz".format(fastq_dir=FASTQ_DIR),
        insertseq = "{combined_insertseq_dir}/insertseqs.combined.fasta".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        i1 = "{combined_insertseq_dir}/insertseqs.combined.fasta.1.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        i2 = "{combined_insertseq_dir}/insertseqs.combined.fasta.2.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        i3 = "{combined_insertseq_dir}/insertseqs.combined.fasta.3.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        i4 = "{combined_insertseq_dir}/insertseqs.combined.fasta.4.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        i5 = "{combined_insertseq_dir}/insertseqs.combined.fasta.rev.1.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR),
        i6 = "{combined_insertseq_dir}/insertseqs.combined.fasta.rev.2.bt2".format(combined_insertseq_dir=COMBINED_INSERTSEQ_DIR)
    output:
        sam = "{insertseq_sam_dir}/{{sample}}.{{pair}}.insertseqs.sam".format(insertseq_sam_dir=INSERTSEQ_SAM_DIR)
    threads:
        config['bowtie2_threads']
    shell:
        "bowtie2 -x {input.insertseq} -U {input.fastq} -S {output.sam} -p {threads} --local --quiet --reorder --no-unal --all"


rule split_genome_alignments:
    input:
        "{genome_sam_dir}/{{sample}}.{{pair}}.genomes.sam".format(genome_sam_dir=GENOME_SAM_DIR)
    output:
        expand("%s/{{sample}}.{{pair}}.{genome}.bam" % SPLIT_GENOME_SAM_DIR, genome=GENOMES)
    shell:
        "python scripts/split_genome_sam.py {input} {output}"


rule split_insertseq_alignments:
    input:
        "{insertseq_sam_dir}/{{sample}}.{{pair}}.insertseqs.sam".format(insertseq_sam_dir=INSERTSEQ_SAM_DIR)
    output:
        expand("%s/{{sample}}.{{pair}}.{insertseq}.bam" % SPLIT_INSERTSEQ_SAM_DIR, insertseq=INSERTSEQS)
    shell:
        "python scripts/split_genome_sam.py {input} {output}"


rule split_sam_by_class:
    input:
        fastq = "{fastq_dir}/{{sample}}.{pair}.fastq.gz".format(fastq_dir=FASTQ_DIR, pair=PAIRS[0]),
        genome_sam_R1 = "{split_genome_sam_dir}/{{sample}}.{pair}.{{genome}}.bam".format(split_genome_sam_dir=SPLIT_GENOME_SAM_DIR, pair=PAIRS[0]),
        genome_sam_R2 = "{split_genome_sam_dir}/{{sample}}.{pair}.{{genome}}.bam".format(split_genome_sam_dir=SPLIT_GENOME_SAM_DIR, pair=PAIRS[1]),
        insertseq_sam_R1 = "{split_insertseq_sam_dir}/{{sample}}.{pair}.{{insertseq}}.bam".format(split_insertseq_sam_dir=SPLIT_INSERTSEQ_SAM_DIR, pair=PAIRS[0]),
        insertseq_sam_R2 = "{split_insertseq_sam_dir}/{{sample}}.{pair}.{{insertseq}}.bam".format(split_insertseq_sam_dir=SPLIT_INSERTSEQ_SAM_DIR, pair=PAIRS[1]),
        class_R1 = "{class_dir}/{{sample}}.{pair}.class.gz".format(class_dir=CLASS_DIR, pair=PAIRS[0]),
        class_R2 = "{class_dir}/{{sample}}.{pair}.class.gz".format(class_dir=CLASS_DIR, pair=PAIRS[1])

    output:
        genome_sam_stats = "{class_genome_sam_dir}/{{sample}}.{{genome}}.{{insertseq}}.stats".format(class_genome_sam_dir=CLASS_GENOME_SAM_DIR),
        insert_sam_stats = "{class_insertseq_sam_dir}/{{sample}}.{{genome}}.{{insertseq}}.stats".format(class_insertseq_sam_dir=CLASS_INSERTSEQ_SAM_DIR)

    shell:
        "python scripts/split_sam_by_class.py "
        "{input.fastq} "
        "{input.genome_sam_R1} {input.genome_sam_R2} "
        "{input.insertseq_sam_R1} {input.insertseq_sam_R2} "
        "{input.class_R1} {input.class_R2} "
        "{output.genome_sam_stats} {output.insert_sam_stats}"


rule sam_merge_taxonomy:
    input:
        taxonomy_db_dir = config['taxonomy_db_dir'],
        genome_sam_stats = "{class_genome_sam_dir}/{{sample}}.{{genome}}.{{insertseq}}.stats".format(class_genome_sam_dir=CLASS_GENOME_SAM_DIR),
    output:
        "{sam_merged_taxonomy_dir}/{{sample}}.{{genome}}.{{insertseq}}.bam".format(sam_merged_taxonomy_dir=SAM_MERGED_TAXONOMY)
    shell:
        "python scripts/sam_merge_taxonomy.py {input.genome_sam_stats} {output} {input.taxonomy_db_dir}"


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
    params:
        depth_cutoff = 5
    shell:
        "python scripts/call_peaks.py {input.bam} {output} {params.depth_cutoff}"


rule plot_peaks:
    input:
        "{peaks_read_depth_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peaks_read_depth_dir=PEAKS_READ_DEPTHS_DIR)
    output:
        "{peak_plots_dir}/{{sample}}.{{genome}}.{{insertseq}}".format(peak_plots_dir=PEAKS_PLOTS_DIR)
    shell:
        "Rscript scripts/create_peak_graphs.R {input} {output}"

