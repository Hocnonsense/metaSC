# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-06 21:57:58
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-10 12:56:59
 * @FilePath: /metaSC/PyLib/reader/read_outputs.py
 * @Description:
        checkM, gtdbtk, iRep, contig_depths, fasta
"""

import os
from typing import Any, Callable, Dict, List, TextIO, Tuple

from Bio import SeqIO
from numpy import nan

from PyLib.PyLibTool.file_info import verbose_import


verbose_import(__name__, __doc__)


def read_nan(dtype: type, numstr: str, *nanstr) -> Any:
    """read args, if arg is nan, return nan
    * @param dtype     type of return
    * @param numstr    number to transfer
                       if numstr in nan: num = str
    * @return type, nan
    """
    if numstr in nanstr:
        return nan
    return dtype(numstr)


def iRep(text: TextIO) -> list:
    """read -o.tsv
     TODO: Sometimes read several .sam. How to solve it?
    @return:
         binId | index of replication (iRep) | un-filtered index of replication (iRep) | raw index of replication (no GC bias correction) | r^2 | coverage | % windows passing filter | fragments/Mbp | GC bias | GC r^2
         n/a will be recorded to np.nan
    """

    def read_values(
        text: TextIO,
        head: str,
        value_type: type,
        recode_func: Callable[[str, Any], None],
    ):
        """
            ../03_modify/7_final//metabat2_90_90.5.fa	n/a
            #
        --> ## r^2
            # genome	E001.sam
            ../03_modify/7_final//17.fa	n/a
            ../03_modify/7_final//4.fa	0.907207353884764
            ../03_modify/7_final//55_sub22.fa	0.9151732818669122
        """
        line = text.readline()
        assert line[: len(head) + 3] + "## " + head
        line = text.readline()  # genome name
        for line in text:
            if line[0] == "#":
                return line
            genome_dir, raw_value = line.strip().split()
            # ../03_modify/7_final//17.fa	n/a
            binId = os.path.basename(os.path.splitext(genome_dir)[0])
            value = read_nan(value_type, raw_value, "n/a", "False")
            recode_func(binId, value)

    # collect data to this map
    irmap = {}
    irvalues = {
        "index of replication": float,
        "un-filtered index of replication": float,
        "raw index of replication": float,
        "r^2": float,
        "coverage": float,
        "% windows passing filter": float,
        "fragments/Mbp": int,
        "GC bias": float,
        "GC r^2": float,
    }
    heads_num = len(irvalues)

    def index_lambda(binId: str, value: Any, index: int) -> None:
        if binId not in irmap:
            irmap[binId] = [nan for _ in range(heads_num)]
        irmap[binId][index] = value

    for index, (name, dtype) in enumerate(irvalues.items()):
        read_values(text, name, dtype, lambda p0, p1: index_lambda(p0, p1, index))
    return [(binId, *values) for binId, values in irmap.items()]


def check_head(head: str):
    """If the {@param looks like a line, then we confirm it is a .depth file
       from jgi_summarize_bam_contig_depths}.
    * @return: [(index, bam_file_name) for bam files]
               if is not a depth file, return None
    """
    head_list = head.strip().split()
    if (
        head_list[0] == "contigName"
        and head_list[1] == "contigLen"
        and head_list[2] == "totalAvgDepth"
        and head_list[3] == head_list[4][:-4]
        and len(head_list) % 2 == 1
    ):
        return [(i, head_list[i]) for i in range(3, len(head_list), 2)]


class ctgDepthSample:
    """
    * @description: help to get depth or given depth of given sample
    * @useage: ctgDepthSample((sample_list, ctg_depth), "totalAvgDepth")[contigName]
    """

    def __init__(self, contig_depths, sample: str, isvar=False) -> None:
        sample_list, ctg_depth = contig_depths
        self.ctg_depth = ctg_depth
        if sample == "length":
            self.i = 0
            self.j = 0
        elif sample == "totalAvgDepth":
            self.i = 0
            self.j = 1
        else:
            self.i = 2 if isvar else 1
            self.j = sample_list.index[sample]

    def __getitem__(self, contigName):
        return self.ctg_depth[contigName][self.i][self.j]


def jgi_depths(
    text: TextIO,
) -> Tuple[List[str], Dict[str, Tuple[Tuple[int, float], List[float], List[float]]]]:
    """Read .depth generated from jgi_summarize_bam_contig_depths
    * @return (
           sample_list: [sample names, ]
           ctg_depth: {contigName: (
               (length, totalAvgDepth),
               [depth in each sample, ], [depth-var in each sample, ])})
    """
    # DEPTH_META = [(0, "contigName"), (1, "contigLen"), (2, "totalAvgDepth")]
    (sample_list, ctg_depth) = ([], {})
    head = next(text)
    sample_index = check_head(head)
    if sample_index is None:
        raise Exception("file is not a depth file")
    sample_list = [i[1] for i in sample_index]
    # sample_index = DEPTH_META + sample_index
    # print(sample_list)
    for line in text:
        values = line.strip().split()
        ctg_depth[values[0]] = (
            (
                int(float(values[1])),
                float(values[2]),
            ),  # to avoid "base 10: '1.2458e+06'"
            [float(values[i[0]]) for i in sample_index],
            [float(values[i[0] + 1]) for i in sample_index],
        )
    return (sample_list, ctg_depth)


def fasta(text):
    """Read all seqs from fasta file
    * @param {TextIO} text
    * @return {dict} seqs: dict -> {record.id: record.seq}
           *NOTE*: record.seq: Bio.Seq.Seq
    """
    return ((record.id, record.seq) for record in SeqIO.parse(text, "fasta"))
