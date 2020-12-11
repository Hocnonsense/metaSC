# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-06 21:57:58
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-11 22:31:53
 * @FilePath: /HScripts/Python/mylib/biotool/read_outputs.py
 * @Description:
        checkM, gtdbtk, iRep, contig_depths, fasta
"""

import os
from io import StringIO
from sys import stderr
from typing import Any, Callable

from Bio import SeqIO
from mylib.biotool.checkm.reload import reload_checkMOutput
from numpy import nan

checkM = reload_checkMOutput


def gtdbtk(text: StringIO) -> list:
    """ read --out_dir/gtdbtk.bac120.summary.tsv
                     ./gtdbtk.ar122.summary.tsv
       @return:
            user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, aa_percent, translation_table, red_value, warnings
       @warnings: some line are need to deal with.
    """
    print(__doc__, file=stderr)
    gtdbtklist = []
    #titles = text.readline().strip().split("\t")
    text.readline()
    print("""
    raise FutureWarning("Some lines are still raw. please check carefully")
    #""", file=stderr)
    for line in text:
        user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, aa_percent, translation_table, red_value, warnings = line[:-1].split("\t")
        gtdbtklist.append([
            user_genome, classification, fastani_reference, fastani_reference_radius, fastani_taxonomy, fastani_ani, fastani_af, closest_placement_reference, closest_placement_taxonomy, closest_placement_ani, closest_placement_af, pplacer_taxonomy, classification_method, note, other_related_references, aa_percent, translation_table, red_value, warnings
        ])
    return gtdbtklist


def read_nan(dtype: type, numstr: str, *nanstr) -> Any:
    """ read args, if arg is nan, return nan
     * @param dtype     type of return
     * @param numstr    number to transfer
                        if numstr in nan: num = str
     * @return type, nan
    """
    if numstr in nanstr:
        return nan
    return dtype(numstr)


def iRep(text: StringIO) -> list:
    """ read -o.tsv
        TODO: Sometimes read several .sam. How to solve it?
       @return:
            binId | index of replication (iRep) | un-filtered index of replication (iRep) | raw index of replication (no GC bias correction) | r^2 | coverage | % windows passing filter | fragments/Mbp | GC bias | GC r^2
            n/a will be recorded to np.nan
    """
    print(__doc__, file=stderr)

    def read_values(text: StringIO, head: str, value_type: type,
                    recode_func: Callable[[str, Any], None]) -> str:
        '''
            ../03_modify/7_final//metabat2_90_90.5.fa	n/a
            #
        --> ## r^2
            # genome	E001.sam
            ../03_modify/7_final//17.fa	n/a
            ../03_modify/7_final//4.fa	0.907207353884764
            ../03_modify/7_final//55_sub22.fa	0.9151732818669122
        '''
        line = text.readline()
        assert line[:len(head) + 3] + "## " + head
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
        'index of replication'            : float,
        'un-filtered index of replication': float,
        'raw index of replication'        : float,
        'r^2'                             : float,
        'coverage'                        : float,
        '% windows passing filter'        : float,
        'fragments/Mbp'                   : int  ,
        'GC bias'                         : float,
        'GC r^2'                          : float,
    }
    heads_num = len(irvalues)

    def index_lambda(binId: str, value: Any, index: int) -> None:
        if binId not in irmap:
            irmap[binId] = [nan for _ in range(heads_num)]
        irmap[binId][index] = value

    for index, (name, dtype) in enumerate(irvalues.items()):
        read_values(text, name, dtype, lambda p0, p1: index_lambda(p0, p1, index))
    return [[binId] + values for binId, values in irmap.items()]


def check_head(head: str) -> bool:
    """ If the {@param looks like a line, then we confirm it is a .depth file
        from jgi_summarize_bam_contig_depths}.
     * @return: [(index, bam_file_name) for bam files]
                if is not a depth file, return None
    """
    head_list = head.strip().split()
    if \
            head_list[0] == "contigName" and \
            head_list[1] == "contigLen" and \
            head_list[2] == "totalAvgDepth" and \
            head_list[3] == head_list[4][:-4] and \
            len(head_list) % 2 == 1:
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


def contig_depths(text: StringIO) -> dict:
    """ Read .depth generated from jgi_summarize_bam_contig_depths
     * @return (
            sample_list: list -> [sample names, ]
            ctg_depth: dict -> {
                contigName: (
                    (length, totalAvgDepth),
                    [depth in each sample, ],
                    [depth-var in each sample]
                )
            }
        )
    """
    print(__doc__, file=stderr)
    # DEPTH_META = [(0, "contigName"), (1, "contigLen"), (2, "totalAvgDepth")]
    (sample_list, ctg_depth) = ([], {})
    head = text.readline()
    sample_index = check_head(head)
    sample_list = [i[1] for i in sample_index]
    #sample_index = DEPTH_META + sample_index
    #print(sample_list)
    for line in text:
        values = line.strip().split()
        ctg_depth[values[0]] = (
            (int(float(values[1])), float(values[2])),
            [float(values[i[0]]) for i in sample_index],
            [float(values[i[0] + 1]) for i in sample_index]
        )
    return (sample_list, ctg_depth)


def fasta(text: StringIO) -> dict:
    """ Read all seqs from fasta file
     * @param {StringIO} text
     * @return {dict} seqs: dict -> {record.id: record.seq}
            *NOTE*: record.seq: Bio.Seq.Seq
    """
    print(__doc__, file=stderr)
    seqs = {}
    for record in SeqIO.parse(text, "fasta"):
        seqs[record.id] = record.seq
    return seqs
