#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   fasta-kmer-split.py
@Time    :   2021/03/05 09:44:27
@Author  :   Bai Shenglong 
@Version :   1.0
@Contact :   slbai01@foxmail.com
@License :   MIT
@Desc    :   None
'''

import argparse, textwrap
from Bio import SeqIO
parser = argparse.ArgumentParser(description='')
parser.add_argument("--input_genome", metavar='\b', type=str, required=True, help=textwrap.dedent('''\
	input genome file
	'''))
parser.add_argument("--input_kmer", metavar='', type=int, required=False, default=23, help="k-mer size, default: 23")
parser.add_argument("-v", "--verbosity", action="count", default=0)
args = parser.parse_args()

dir_kmer_times = {}
dir_kmer_info = {}
def process_one_fasta(fasta_id, fasta_seq, kmer_size):
	for x in range(0, len(fasta_seq) - kmer_size):
		kmer_seq = fasta_seq[x : x + kmer_size]
		kmer_start = x
		try:
			dir_kmer_times[kmer_seq] = dir_kmer_times[kmer_seq] + 1
		except:
			dir_kmer_times[kmer_seq] = 1
			dir_kmer_info[kmer_seq] = "{}\t{}\t{}".format(kmer_seq, fasta_id, kmer_start)

kmer_size = args.input_kmer
for record in SeqIO.parse(args.input_genome, "fasta"):
	fasta_id = record.id
	fasta_seq = record.seq.upper()
	fasta_rev_comp_seq = fasta_seq.reverse_complement()
	process_one_fasta(fasta_id, fasta_seq, kmer_size)
	process_one_fasta(fasta_id, fasta_rev_comp_seq, kmer_size)

for kmer_seq, kmer_times in dir_kmer_times.items():
	if kmer_times == 1:
		print(dir_kmer_info[kmer_seq])

