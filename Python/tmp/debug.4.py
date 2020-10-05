# -*- coding: utf-8 -*-
"""
    @Editor: Lv Yongxin
    @Date: 2020/08/27 21:16:31
    @LastEditors: Hwrn
    @LastEditTime: 2020-09-15 20:29:11
    @FilePath: /mylib/debug.py
    @Description: 对prodigal或者refinem的结果进行修剪，小于33bp的序列都删掉了
    @Version: 1.0.1
"""


from sys import argv
import sys
from Bio import SeqIO


__help = """
    trim results from prodigal. sequence less then 33 bp are discarded.
"""

if len(argv != 3):
    print(__help)
    exit(1)

sc, prefix, suffix = argv

# get gene length and deside which to discard
dict_gene_needed = {}
try:
    with open(prefix+'.faa') as fin, \
            open(suffix+'.faa', 'w') as fout:
        for record in SeqIO.parse(fin, 'fasta'):
            des = str(record.description)
            sid = str(record.id)
            seq = str(record.seq)
            if len(seq) >= 33:
                dict_gene_needed[sid] = 1
                fout.write('>'+des+'\n'+seq+'\n')
except IOError:
    print("must provide amino acid seq file")
    sys.exit()

try:
    with open(prefix+'.fna') as fin, open(suffix+'.fna', 'w') as fout:
        for record in SeqIO.parse(fin, 'fasta'):
            des = str(record.description)
            sid = str(record.id)
            seq = str(record.seq)
            if sid in dict_gene_needed:
                dict_gene_needed[sid] = 1
                fout.write('>'+des+'\n'+seq+'\n')
except IOError:
    print("nucleotide seq file not found,pass")

try:
    with open(prefix+'.gff') as fin, open(suffix+'.gff', 'w') as fout:
        for line in fin:
            if line[0] != '#':
                temp = line.strip().split('\t')
                scaffold = temp[0]
                s_count = temp[8].split(';')[0].split('_')[1]
                sid = scaffold+'_'+s_count
                if sid in dict_gene_needed:
                    fout.write(line)
except IOError:
    print("nucleotide seq file not found,pass")
