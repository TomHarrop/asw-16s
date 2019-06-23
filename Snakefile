#!/usr/bin/env python3

import multiprocessing

#############
# FUNCTIONS #
#############

def indiv_to_fastq(wildcards):
    return(
        {'r1': indiv_to_r1[wildcards.indiv],
         'r2': indiv_to_r2[wildcards.indiv]})

###########
# GLOBALS #
###########

bbmap_container = 'shub://TomHarrop/singularity-containers:bbmap_38.50b'
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.8beta'
bracken_container = 'shub://TomHarrop/singularity-containers:bracken_2.2'
bioc_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
biopython_container= 'shub://TomHarrop/singularity-containers:biopython_1.73'

r1_path = 'data/raw_full/4826-{indiv}-0-1_S{sindiv}_L001_R1_001.fastq.gz'
r2_path = 'data/raw_full/4826-{indiv}-0-1_S{sindiv}_L001_R2_001.fastq.gz'

########
# MAIN #
########

all_indivs = sorted(set(glob_wildcards(r1_path).indiv))

# exclude negatives from kraken
neg_indivs = ['211',
              '212',
              '213',
              '214',
              '215',
              '216',
              '217',
              '218',
              '219']

kraken_indivs = [x for x in all_indivs
                 if x not in neg_indivs]

indiv_to_r1 = {
    x: r1_path.format(indiv=x, sindiv=int(x))
    for x in all_indivs}
indiv_to_r2 = {
    x: r2_path.format(indiv=x, sindiv=int(x))
    for x in all_indivs}

#########
# RULES #
#########

rule target:
    input:
        'output/010_trimmed/merge_stats.txt',
        expand('output/010_trimmed/{indiv}/readlength.txt',
               indiv=all_indivs),
        'output/040_phyloseq/filtered_reads_per_indiv.pdf',
        'output/040_phyloseq/alpha_diversity.pdf',
        'output/040_phyloseq/ordination.pdf'


rule ordination:
    input:
        phyloseq = 'output/040_phyloseq/ps_filtered.Rds',
    output:
        plot = 'output/040_phyloseq/ordination.pdf'
    log:
        'output/logs/040_phyloseq/ordination.log'
    singularity:
        bioc_container
    script:
        'src/ordination.R'


rule alpha_diversity:
    input:
        phyloseq = 'output/040_phyloseq/ps_filtered.Rds',
    output:
        plot = 'output/040_phyloseq/alpha_diversity.pdf'
    log:
        'output/logs/040_phyloseq/alpha_diversity.log'
    singularity:
        bioc_container
    script:
        'src/alpha_diversity.R'


rule filtered_reads_per_indiv:
    input:
        phyloseq = 'output/040_phyloseq/ps_filtered.Rds',
    output:
        plot = 'output/040_phyloseq/filtered_reads_per_indiv.pdf'
    log:
        'output/logs/040_phyloseq/filtered_reads_per_indiv.log'
    singularity:
        bioc_container
    script:
        'src/filtered_reads_per_indiv.R'


rule filter_otus:
    input:
        phyloseq = 'output/040_phyloseq/ps.Rds',
    output:
        phyloseq = 'output/040_phyloseq/ps_filtered.Rds',
    log:
        'output/logs/040_phyloseq/filter_otus.log'
    singularity:
        bioc_container
    script:
        'src/filter_otus.R'


rule construct_phyloseq:
    input:
        mpa_files = expand(
            'output/030_bracken/{indiv}/bracken_report_mpa.txt',
            indiv=kraken_indivs),
        sample_catalog = 'data/sample_catalog.csv',
        indiv_mapping = 'data/ogbf_sample_info.csv'
    output:
        phyloseq = 'output/040_phyloseq/ps.Rds',
        parsed_data = 'output/040_phyloseq/parsed_data.csv'
    log:
        'output/logs/040_phyloseq/construct_phyloseq.log'
    singularity:
        bioc_container
    script:
        'src/construct_phyloseq.R'


rule mpa_report:
    input:
        'output/020_kraken/{indiv}/kraken_report_bracken.txt'
    output:
        'output/030_bracken/{indiv}/bracken_report_mpa.txt'
    log:
        'output/logs/030_bracken/{indiv}_bracken-mpa.log'
    singularity:
        bracken_container
    shell:
        'kreport2mpa.py '
        '-r {input} '
        '-o {output} '
        '&> {log}'

rule bracken:
    input:
        report = 'output/020_kraken/{indiv}/kraken_report.txt',
        db = 'data/20190614-silva',
    params:
        readlength = '463',
        threshold = '1',
        level = 'G'
    output:
        report = 'output/030_bracken/{indiv}/bracken_report.txt',
        crap_file = 'output/020_kraken/{indiv}/kraken_report_bracken.txt'
    log:
        'output/logs/030_bracken/{indiv}_bracken.log'
    singularity:
        bracken_container
    shell:
        'bracken '
        '-d {input.db} '
        '-i {input.report} '
        '-o {output.report} '
        '-r {params.readlength} '
        '-t {params.threshold} '
        '-l {params.level} '
        '&> {log} '

rule kraken:
    input:
        fq = 'output/010_trimmed/{indiv}/merged.fastq.gz',
        db = 'data/20190614-silva'
    output:
        out = 'output/020_kraken/{indiv}/kraken_out.txt',
        report = 'output/020_kraken/{indiv}/kraken_report.txt'
    log:
        'output/logs/020_kraken/{indiv}_kraken.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--report-zero-counts '
        # '--use-mpa-style '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.fq} '
        '&> {log}'

rule merged_stats:
    input:
        fq = expand('output/010_trimmed/{indiv}/merged.fastq.gz',
                    indiv=all_indivs)
    output:
        'output/010_trimmed/merge_stats.txt'
    params:
        inline = lambda wildcards,input: ','.join(input.fq)
    log:
        'output/logs/010_trimmed/statswrapper.log'
    singularity:
        bbmap_container
    shell:
        'statswrapper.sh '
        'in={params.inline} '
        '> {output} '
        '2> {log}'

rule readlength:
    input:
        'output/010_trimmed/{indiv}/merged.fastq.gz'
    output:
        'output/010_trimmed/{indiv}/readlength.txt'
    log:
        'output/logs/010_trimmed/{indiv}_readlength.log'
    singularity:
        bbmap_container
    shell:
        'readlength.sh '
        'in={input} '
        'out={output} '
        'bin=1 '
        'max=1000 '
        'nzo=f '
        '2> {log}'


rule trim_merge:
    input:
        unpack(indiv_to_fastq),
        # barcodes = 'output/barcodes.fasta' # not currently working
    output:
        fq = 'output/010_trimmed/{indiv}/merged.fastq.gz',
        stats = 'output/010_trimmed/{indiv}/trimstats.txt',
        ihist = 'output/010_trimmed/{indiv}/ihist.txt',
    params:
        ref = '/phix174_ill.ref.fa.gz',
        adaptors = '/adapters.fa'
    log:
        bbduk1 = 'output/logs/010_trimmed/{indiv}_bbduk1.log',
        bbduk2 = 'output/logs/010_trimmed/{indiv}_bbduk2.log',
        bbmerge = 'output/logs/010_trimmed/{indiv}_bbmerge.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbmap_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        '-Xmx12g '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'stats={output.stats} '
        'forcetrimmod=5 '
        'k=27 '
        'mink=10 '
        'ktrim=l tpe tbo '
        'editdistance=2 '
        'editdistance2=2 '
        # 'barcodefilter=t '
        # 'barcodes={input.barcodes} '
        'ref={params.ref} '
        '2> {log.bbduk1} '
        ' | '
        'bbmerge.sh '
        '-Xmx12g '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'maxlength=500 '
        'mininsert=400 '
        'adapter={params.adaptors} '
        'ihist={output.ihist} '
        'strict=t '
        '2> {log.bbmerge}'
        ' | '
        'bbduk.sh '
        'threads={threads} '
        '-Xmx12g '
        'in=stdin.fastq '
        'int=t '
        'maxns=0 '
        'literal='
        'AAAAAAAAAA,'
        'CCCCCCCCCC,'
        'GGGGGGGGGG,'
        'TTTTTTTTTT '
        'maskmiddle=f '
        'out={output.fq} '
        'ziplevel=9 '
        '2> {log.bbduk2} '

rule generate_bc_fasta:
    input:
        barcodes = 'data/ogbf_sample_info.csv'
    output:
        barcodes = 'output/barcodes.fasta'
    log:
        'output/logs/generate_bc_fasta.log'
    singularity:
        biopython_container
    script:
        'src/generate_bc_fasta.py'

