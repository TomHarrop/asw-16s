#!/usr/bin/env python3

import csv
import logging
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

# set up log
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename=snakemake.log[0],
    level=logging.INFO)

bc_file = snakemake.input['barcodes']
bc_fasta = snakemake.output['barcodes']

# dev
# bc_file ='data/ogbf_sample_info.csv'
# bc_fasta = 'test.fasta'

seq_list = []

with open(bc_file, 'rt') as csv_file:
    csv_reader = csv.reader(csv_file)
    header = next(csv_reader)
    # for x in csv_reader:
    #     if x[0] == '4826-159':
    #         row = x    
    for row in csv_reader:
        my_seq = f'{Seq.Seq(row[3]).reverse_complement()}+{row[5]}'
        my_req = SeqRecord.SeqRecord(Seq.Seq(my_seq),
                                     id=row[0],
                                    description='')
        seq_list.append(my_req)

SeqIO.write(seq_list, bc_fasta, 'fasta')
