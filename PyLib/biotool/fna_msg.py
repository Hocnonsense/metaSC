# -*- coding: utf-8 -*-
"""
 * @Date: 2020-12-11 10:22:23
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 12:01:10
 * @FilePath: /metaSC/PyLib/biotool/fna_msg.py
 * @Description:
        Get message from a fna file.
"""

from typing import Iterable, List, Tuple, Dict

# from PyLib.reader.read_outputs import jgi_depths
from PyLib.PyLibTool.file_info import verbose_import


verbose_import(__name__, __doc__)


def seq_GC_len(seqs: Iterable[Tuple[str, str]]) -> Dict[str, Tuple[int, int, int]]:
    """Read fasta files.
    * @param {Dict} seqs: Iterable[Tuple[str, str]] -> {record.id: record.seq}
    * @return {dict} {ctg_name: [GC count, genome size]}
    """
    GC_len = {}
    for id, seq in seqs:
        seq = seq.upper().replace("-", "").replace(" ", "").replace("\n", "")
        gc_count = seq.count("G") + seq.count("C")
        seq_len = len(seq)
        N_count = seq_len - gc_count - seq.count("A") - seq.count("T")
        GC_len[id] = gc_count, seq_len, N_count
    return GC_len


def length_NL(seqs: Iterable[Tuple[str, str]], pcg=0.5) -> Tuple[int, int]:
    """Calculate N50, L50 or N{pcg}, L{pcg} for given pcg (total is 100)
    * @param {dict} seqs: Iterable[Tuple[str, str]] -> {record.id: record.seq}
    * @param {int} pcg: threashold of total length percentage (100)
    * @return {(int, int)} (N50, L50)
    """
    seqs_len = [len(seq) for seq in seqs]
    seqs_len.sort(reverse=True)
    pcg_len = sum(seqs_len) * pcg
    for i, length in enumerate(seqs_len):
        pcg_len -= length
        if pcg_len <= 0:
            break
    return length, i


def statistic_fna(seqs: Iterable[Tuple[str, str]]) -> Tuple:
    """{
        'header': ['SeqNumbers', 'MaxLength', 'GenomeSize',
                   'GC', 'N50', 'L50']
    }"""
    GC_len = sorted(seq_GC_len(seqs).values(), reverse=True)
    SeqNumbers = len(GC_len)
    MaxLength = GC_len[0][1]

    total_len = sum(i[1] for i in GC_len)
    GenomeSize = total_len
    GC = sum(i[0] for i in GC_len) / sum(i[1] - i[2] for i in GC_len)

    pcg_len = total_len * 50 / 100
    for i, (gc, length, _) in enumerate(GC_len):
        pcg_len -= length
        if pcg_len <= 0:
            N50 = length
            L50 = i + 1  # start from 1
            break

    return SeqNumbers, MaxLength, GenomeSize, GC, N50, L50


def seq_total_depth(
    ctg_depth: dict,
    seqs: Iterable[Tuple[str, str]],
) -> Tuple[float, List[float]]:
    """
    * @description:
       It is difficult to tell how depth is calculated by *jgi_summarize_bam_contig_depths*
       https://bitbucket.org/berkeleylab/metabat/issues/48/jgi_summarize_bam_contig_depths-coverage
           #L699: int32_t ignoreEdges = 0

           #L715: end = contigLength
           #L717: start = ignoreEdges;
           #L718: end = end - ignoreEdges;

           #L721: int32_t adjustedContigLength = end - start

           #L728: avgDepth += depthCounts.baseCounts[i] * (hasWeights ? weights[i] : 1.0);
           #L732: avgDepth = avgDepth / adjustedContigLength;
           total_bases = sum(avgDepth * length)
           total_depth = total_bases / total_length
    * @param {Iterable} seqs: [scaffold_name, ]
    * @param {dict} ctg_depth: {
           contigName: (
               (length, totalAvgDepth),
               [depth in each sample, ], [depth-var in each sample, ]
           )
       } from jgi_depths
    * @return
       (totalAvgDepth, [depth in each sample]) of given seqs
    """
    (totalLength, totalBases) = (0, 0.0)
    for values in ctg_depth:
        sample_len = len(ctg_depth[values][1])
        break  # get the length and leave
    sampleBases = [0.0] * sample_len
    for contigName, _ in seqs:
        (seq_length, totalAvgDepth), sample_depth, sample_depth_var = ctg_depth[
            contigName
        ]
        totalLength += seq_length
        totalBases += seq_length * totalAvgDepth
        for i, depth in enumerate(sample_depth):
            sampleBases[i] += depth * seq_length
    return (totalBases / totalLength, [bases / totalLength for bases in sampleBases])
