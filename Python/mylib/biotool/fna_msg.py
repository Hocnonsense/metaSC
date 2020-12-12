# -*- coding: utf-8 -*-
"""
 * @Date: 2020-12-11 10:22:23
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-11 17:47:42
 * @FilePath: /HScripts/Python/mylib/biotool/fna_msg.py
 * @Description:
        Get message from a fna file.
"""

from collections import OrderedDict
from sys import stderr
from typing import Iterable, Tuple


def seq_GC_len(seqs: dict) -> dict:
    """ Read fasta files.
     * @param {dict} seqs: dict -> {record.id: record.seq}
     * @return {dict} {ctg_name: [GC count, genome size]}
    """
    print(__doc__, file=stderr)
    GC_len = {}
    for id, seq in seqs.items():
        gc_count = seq.count("G") + seq.count("C")
        seq_len = len(seq)
        GC_len[id] = gc_count, seq_len
    return GC_len


def length_NL(seqs: dict, pcg=50) -> Tuple[int, int]:
    """ Calculate N50, L50 or N{pcg}, L{pcg} for given pcg (total is 100)
     * @param {dict} seqs: dict -> {record.id: record.seq}
     * @param {int} pcg: threashold of total length percentage (100)
     * @return {(int, int)} (N50, L50)
    """
    seqs_len = [len(seq) for seq in seqs.items()]
    seqs_len.sort(reverse=True)
    pcg_len = sum(seqs_len) * pcg / 100
    for i, length in enumerate(seqs_len):
        pcg_len -= length
        if pcg_len <= 0:
            return (length, i)


def statistic_fna(seqs: dict) -> dict:
    """ Sequence message: SeqNumbers, MaxLength, GenomeSize, GC, N50, L50
     * @param {dict} seqs: dict -> {record.id: record.seq}
     * @return {dict} fna_msg
    """
    fna_msg = OrderedDict()
    GC_len = sorted(seq_GC_len(seqs).values(), reverse=True)
    fna_msg["SeqNumbers"] = len(GC_len)
    fna_msg["MaxLength"] = GC_len[0][1]

    total_len = sum(i[1] for i in GC_len)
    fna_msg["GenomeSize"] = total_len
    fna_msg["GC"] = sum(i[0] for i in GC_len) / total_len

    pcg_len = total_len * 50 / 100
    for i, (gc, length) in enumerate(GC_len):
        pcg_len -= length
        if pcg_len <= 0:
            fna_msg["N50"] = length
            fna_msg["L50"] = i
            break

    return fna_msg


def seq_depth(seqs: Iterable, ctg_depth: dict) -> dict:
    """ Get depth of given seqs subset.
    * @param {Iterable} seqs: [scaffold_name, ]
    * @param {dict} ctg_depth: dict -> {
            contigName: (
                (length, totalAvgDepth),
                [depth in each sample, ],
                [depth-var in each sample]
            )
        } from contig_depths
    * @return {dict} sub_ctg_depth: subset of ctg_depth
    """
    return {contigName: ctg_depth[contigName] for contigName in seqs}


def seq_total_depth(seqs: Iterable, ctg_depth: dict):
    """ Total bases and depth of given seqs
            total_bases = sum(avgDepth * length)
            total_depth = total_bases / total_length
    * @param {Iterable} seqs: [scaffold_name, ]
    * @param {dict} ctg_depth: dict -> {
            contigName: (
                (length, totalAvgDepth),
                [depth in each sample, ],
                [depth-var in each sample]
            )
        } from contig_depths
     * @return {(float, list)} (totalAvgDepth, [depth in each sample]) of given seqs
    """
    (totalLength, totalBases) = (0, 0.0)
    for values in ctg_depth:
        sample_len = len(values[1])
        break  # get the length and leave
    sampleBases = [0.0] * sample_len
    for contigName in seqs:
        (seq_length, totalAvgDepth), sample_depth, sample_depth_var = ctg_depth[contigName]
        totalLength += seq_length
        totalBases += seq_length * totalAvgDepth
        for i, depth in enumerate(sample_depth):
            sampleBases[i] += depth * seq_length
    return (totalBases / totalLength, [bases / totalLength for bases in sampleBases])
