# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-24 10:24:10
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-05 11:49:33
 * @FilePath: /HScripts/Python/seqPipe/x01_length_cut.py
 * @Description:
    01_length_cut.py <in_file.fa> <out_file.fa> <threshold> [-h] [--help]
        @description:   Remove sequence shorter than given threshold.
                        This file is designed to control quality of
                        scaffolds.fa generated by assembly program like
                        megahit, spades, or idba.
        @param:         <in_file.fa>:   input file in FASTA format
                        <out_file.fa>:  output file
                        <threshold>:    threshold, typtically, 500
                                        Any sequence < threshold will be
                                        discard
"""

import os
from sys import argv, stderr
from Bio import SeqIO


def parse_args():
    if "-h" in argv or "--help" in argv or len(argv) == 1:
        print(__doc__, file=stderr)
        exit(0)

    sc, in_file, out_file, number = argv
    in_file = os.path.abspath(os.path.expanduser(in_file))
    out_file = os.path.abspath(os.path.expanduser(out_file))
    num = int(number)
    args = in_file, out_file, num
    print(sc, *args, sep="\n" + " "*4, file=stderr)
    return args


def main():
    in_file, out_file, num = parse_args()

    with open(out_file, "w") as fo \
            , open(in_file) as fi \
            :
        discard_seqs, discard_bases = 0, 0
        for line in SeqIO.parse(fi, 'fasta'):
            if len(line.seq) >= num:
                fo.write('>'+str(line.id)+'\n'+str(line.seq)+'\n')
            else:
                discard_seqs += 1
                discard_bases += len(line.seq)
    print("{seqs_n} reads ({bases_n} bases) are discarded\n".format(
        seqs_n=discard_seqs, bases_n=discard_bases), file=stderr)


if __name__ == "__main__":
    main()
