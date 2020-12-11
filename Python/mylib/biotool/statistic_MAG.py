# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 23:09:57
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-11 17:10:41
 * @FilePath: /HScripts/Python/mylib/biotool/statistic_MAG.py
 * @Description:
    seq number, GC%, genome size from *.fa file
"""

import os
from sys import argv, stderr

from mylib.biotool.read_outputs import checkM, contig_depths, fasta, gtdbtk, iRep
from mylib.biotool.fna_msg import statistic_fna, seq_total_depth


def list_MAGs(MAG_file_path: str, endswith="") -> list:
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


def collect_MAGs_msg(
        MAG_file_path: str,
        depth_file: str = None,
        checkM_file: str = None,
        iRep_file: str = None,
        gtdbtk_file: str = None
):
    """
     * @return {(list, dict)} -> {
         [title for each rows],
         {
             MAG: [values for each rows]
         }
     * }
    """
    msg_title = []
    MAGs_msg = {}

    # check dir
    MAG_file_path = os.path.abspath(os.path.expanduser(MAG_file_path))
    MAGs_list = list_MAGs(MAG_file_path)
    # calculate depth with seqs
    MAG_msg_title = []
    if depth_file:
        with open(depth_file) as text:
            sample_list, ctg_depth = contig_depths(text)

    for MAG in MAGs_list:
        MAG_msg = []
        MAG_file = os.path.join(MAG_file_path, MAG)
        MAG_seqs = fasta(MAG_file)
        fna_msg = statistic_fna(MAG_seqs)
        if not MAG_msg_title:
            MAG_msg_title = list(fna_msg)

        MAGs_msg[MAG] = MAG_msg
        for k in MAG_msg_title:
            MAG_msg.append(fna_msg[k])

        depth_title = []
        if depth_file:
            depth_title = ["totalAvgDepth"] + sample_list
            totalAvgDepth, sampleDepth = seq_total_depth(MAG_seqs, ctg_depth)
            MAG_msg += [totalAvgDepth] + sampleDepth

    msg_title += MAG_msg_title + depth_title

    if iRep_file:
        with open(iRep_file) as text:
            iReplist = iRep(text)
        msg_title += ["index of replication"]
        iRep_MAG_msg = {}
        for msg in iReplist:
            iRep_MAG_msg[msg[0]] = [msg[1]]
        for MAG, MAG_msg in MAGs_msg.items():
            MAG_msg += iRep_MAG_msg.get(MAG, ["n/a"])

    if gtdbtk_file:
        with open(gtdbtk_file) as text:
            gtdbtklist = gtdbtk()
        gtdbtk_title = ["Domain", "Phylum", "Class",
                        "Order", "Family", "Genus", "Species",
                        "ClosestRefGenome", "ANI"]
        gtdbtk_MAG_msg = {}
        for msg in gtdbtklist:
            gtdbtk_MAG_msg[msg[0]] = msg[1].split(";") + msg[2:4]
        for MAG, MAG_msg in MAGs_msg.items():
            MAG_msg += gtdbtk_MAG_msg.get(os.path.basename(MAG),
                                          ["n/a"] * len(gtdbtk_title))
        msg_title += gtdbtk_title

    if checkM_file:
        with open(checkM_file) as text:
            ckmap = checkM(checkM_file)
        msg_title += ["Marker_lineage", "Completeness",
                      "Contamination", "Strain_heterogeneity"]
        checkM_MAG_msg = {}
        for Bin_Id, values in ckmap.items():
            checkM_MAG_msg[Bin_Id] = [values[1]] + values[-3:]
        for MAG, MAG_msg in MAGs_msg.items():
            MAG_msg += checkM_MAG_msg.get(os.path.basename(MAG),
                                          ["n/a"] * 4)

    return (msg_title, MAGs_msg)


if __name__ == "__main__":
    args = argv[1:]
    # TODO: add func to use it with bash
