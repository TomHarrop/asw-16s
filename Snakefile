#!/usr/bin/env python3


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

r1_path = 'data/raw/4826-{indiv}-0-1_S{sindiv}_L001_R1_001.fastq.gz'
r2_path = 'data/raw/4826-{indiv}-0-1_S{sindiv}_L001_R2_001.fastq.gz'

########
# MAIN #
########

all_indivs = sorted(set(glob_wildcards(r1_path).indiv))

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
        'output/010_trimmed/merge_stats.txt'

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

rule trim_merge:
    input:
        unpack(indiv_to_fastq)
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
    singularity:
        bbmap_container
    shell:
        'bbduk.sh '
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
        'ref={params.ref} '
        '2> {log.bbduk1} '
        ' | '
        'bbmerge.sh '
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
