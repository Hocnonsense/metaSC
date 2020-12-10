# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 23:09:57
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-10 20:37:27
 * @FilePath: /HScripts/Python/mylib/biotool/statistic_MAG.py
 * @Description:
    seq number, GC%, genome size from *.fa file
"""

from io import StringIO
import os
from sys import stderr
from typing import Tuple
from Bio import SeqIO


def list_MAGs(MAG_file_path: str, endswith = "") -> list:
    """
    * @description: list name of MAGs in given path
    * @param {str} MAG_file_path
    * @param {str} endswith
    * @return {list} MAGs_list: [name of MAGs (with "endswith")]
    """
    print(__doc__, file=stderr)
    MAGs_list = []
    for MAG_file in sorted(os.listdir(MAG_file_path)):
        if MAG_file.endswith(endswith):
            MAGs_list.append(MAG_file)
    return MAGs_list


def get_MAG_ctgs(MAG_file: Tuple[str, StringIO]) -> dict:
    """ now, read MAG's fasta files
     * @return {dict} MAG_dict: [scaffold_name, ] of given MAG
    """
    print(__doc__, file=stderr)
    MAG_ctgs = []
    for record in SeqIO.parse(MAG_file, "fasta"):
        MAG_ctgs.append(record.name)
    return MAG_ctgs


def get_ctg_msg(fasta_file: Tuple[str, StringIO]) -> list:
    """ Read fasta files.
     * @param MAG_file_path: path of MAG file or scaffold.fa or IO.
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


def get_ctg_depth(ctgs: list, ctg_depth: dict) -> list:
    """ Get depth of given MAG.
    * @param {list} ctgs: [scaffold_name, ]
    * @param {dict} ctg_depth: dict -> {
            contigName: (
                (length, totalAvgDepth),
                [depth in each sample, ],
                [depth-var in each sample]
            )
        } from contig_depths
    * @return {dict} sub_ctg_depth: subset of ctg_depth
    """
    return {contigName: ctg_depth[contigName] for contigName in ctgs}


def sum_ctg_depth(ctgs: list, ctg_depth: dict) -> tuple:
    """ Calculate total depth of given MAG (in all MAGs).
    * @param {list} ctgs: [scaffold_name, ]
    * @param {dict} ctg_depth: dict -> {
            contigName: (
                (length, totalAvgDepth),
                [depth in each sample, ],
                [depth-var in each sample]
            )
        } from contig_depths
    * @return {tuple} MAG_depth_sum: (
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
    for contigName in ctgs:
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


def collect_MAG_segs(
        MAG_file_path: str,
        ):
    pass
