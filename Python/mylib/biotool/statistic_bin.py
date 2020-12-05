# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 23:09:57
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-11-14 15:56:19
 * @FilePath: /HScripts/Python/mylib/biotool/statistic_bin.py
 * @Description:
    seq number, GC%, genome size from *.fa file
"""

from io import StringIO
import os
from sys import stderr
from typing import Tuple
from Bio import SeqIO


def list_bins(bin_file_path: str, endswith = "") -> list:
    """
    * @description: list name of bins in given path
    * @param {str} bin_file_path
    * @param {str} endswith
    * @return {list}: [name of bins (with "endswith")]
    """
    print(__doc__, file=stderr)
    bins_list = []
    for bin_file in sorted(os.listdir(bin_file_path)):
        if bin_file.endswith(endswith):
            bins_list.append(bin_file)
    return bins_list


def get_bin_ctgs(bin_file: Tuple[str, StringIO]) -> dict:
    """ now, read bin's fasta files
     * @return bin_dict: [scaffold_name, ] of given bin
    """
    print(__doc__, file=stderr)
    bin_ctgs = []
    for record in SeqIO.parse(bin_file, "fasta"):
        bin_ctgs.append(record.name)
    return bin_ctgs


def get_ctg_msg(fasta_file: Tuple[str, StringIO]) -> list:
    """ Read fasta files.
     * @param bin_file_path: path of bin file or scaffold.fa or IO.
     * @return {dict} {ctg_name: [genome size, GC%]}
    """
    print(__doc__, file=stderr)
    ctgs_msg = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ctg_name = record.name
        seq = record.seq
        gc_count = seq.count("G") + seq.count("C")
        seq_len = len(seq)
        ctgs_msg[ctg_name] = seq_len, gc_count / seq_len
    return ctgs_msg


def get_bin_depth(bin_ctgs: list, ctg_depth: dict) -> list:
    """ Get depth of given bin.
    * @param {list} bin_dict: [scaffold_name, ] of given bin by get_bin_ctgs
    * @param {dict} ctg_depth: dict -> {
            contigName: (
                (length, totalAvgDepth),
                [depth in each sample, ],
                [depth-var in each sample]
            )
        } from contig_depths
    * @return {*}
    """
    return {contigName: ctg_depth[contigName] for contigName in bin_ctgs}


def sum_bin_depth(bin_ctgs: list, ctg_depth: dict) -> tuple:
    """ Calculate total depth of given bin (in all bins).
    * @param {list} bin_dict: [scaffold_name, ] of given bin by get_bin_ctgs
    * @param {dict} ctg_depth: dict -> {
            contigName: (
                (length, totalAvgDepth),
                [depth in each sample, ],
                [depth-var in each sample]
            )
        } from contig_depths
    * @return {tuple} bin_depth_sum: (
            (length, totalAvgDepth),
            [depth in each sample, ],
            [depth-var in each sample]
        )
    """
    (length, totalAvgDepth) = (0, 0.0)
    for values in ctg_depth:
        sample_len = len(values[1])
        sample_depths = [0.0 for _ in values[1]]
        sample_depths_var = [0.0 for _ in values[1]]
        break  # get the length and leave
    for contigName in bin_ctgs:
        values = ctg_depth[contigName]
        length += values[0][0]
        totalAvgDepth += values[0][1]
        for i in range(sample_len):
            sample_depths[i] += values[1][i]
            sample_depths_var[i] += values[2][i]
    return (
        (length, totalAvgDepth),
        sample_depths,
        sample_depths_var
    )
