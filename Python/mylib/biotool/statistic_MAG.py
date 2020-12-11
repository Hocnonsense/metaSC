# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 23:09:57
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-12-11 23:04:56
 * @FilePath: /HScripts/Python/mylib/biotool/statistic_MAG.py
 * @Description:
    seq number, GC%, genome size from *.fa file
"""

import os
from sys import stderr

from mylib.biotool.fna_msg import seq_total_depth, statistic_fna
from mylib.biotool.read_outputs import (checkM, contig_depths, fasta, gtdbtk,
                                        iRep)


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
        MAG_file_path   : str,
        depth_file      : str = "",
        checkM_file     : str = "",
        iRep_file       : str = "",
        gtdbtk_path     : str = ""
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

    try:
        with open(iRep_file) as text:
            iReplist = iRep(text)
        msg_title += ["index of replication"]
        iRep_MAG_msg = {}
        for msg in iReplist:
            iRep_MAG_msg[msg[0]] = [msg[1]]
        for MAG, MAG_msg in MAGs_msg.items():
            MAG_msg += iRep_MAG_msg.get(os.path.basename(
                os.path.splitext(MAG)[0]), ["n/a"])
    except FileNotFoundError:
        pass

    gtdbtk_title = ["Domain", "Phylum", "Class",
                    "Order", "Family", "Genus", "Species",
                    "ClosestRefGenome", "ANI"]
    gtdbtk_MAG_msg = {}
    for markers in ["ar122", "bac120"]:
        gtdbtk_file = os.path.join(
            gtdbtk_path, "gtdbtk.{markers}.markers_summary.tsv".format(markers=markers))
        try:
            with open(gtdbtk_file) as text:
                gtdbtklist = gtdbtk(text)
            for msg in gtdbtklist:
                gtdbtk_MAG_msg[msg[0]] = msg[1].split(";") + [msg[2], msg[3]]
        except FileNotFoundError:
            print(gtdbtk_file, "not found, pass", file=stderr)
    if gtdbtk_MAG_msg:
        for MAG, MAG_msg in MAGs_msg.items():
            MAG_msg += gtdbtk_MAG_msg.get(os.path.basename(
                os.path.splitext(MAG)[0]), ["n/a"] * len(gtdbtk_title))
        msg_title += gtdbtk_title
    else:
        pass

    try:
        with open(checkM_file) as text:
            ckmap = checkM(checkM_file)
        msg_title += ["Marker_lineage", "Completeness",
                      "Contamination", "Strain_heterogeneity"]
        checkM_MAG_msg = {}
        for Bin_Id, values in ckmap.items():
            checkM_MAG_msg[Bin_Id] = [values[1]] + values[-3:]
        for MAG, MAG_msg in MAGs_msg.items():
            MAG_msg += checkM_MAG_msg.get(os.path.basename(
                os.path.splitext(MAG)[0]), ["n/a"] * 4)
    except FileNotFoundError:
        pass

    return (msg_title, MAGs_msg)


if __name__ == "__main__":
    # args = argv[1:]
    # TODO: add func to use it with bash
    # from mylib.biotool.statistic_MAG import collect_MAGs_msg
    pass
    MAG_file_path = "/home/hwrn/Work/2020-12-Yap/04_MAG/00_raw/{sample}/DAS_Tool_DASTool_bins"
    depth_file = "/home/hwrn/Work/2020-12-Yap/02_assembly/depth/{sample}-jgi.depth"
    checkM_file = "/home/hwrn/Work/2020-12-Yap/04_MAG/01_modify/{sample}/01/checkm.out"
    iRep_file = "/home/hwrn/Work/2020-12-Yap/04_MAG/iRep/{sample}.tsv"
    gtdbtk_file = "/home/hwrn/Work/2020-12-Yap/04_MAG/gtdbtk/{sample}/"

    fo = open("/home/hwrn/Work/2020-12-Yap/Analyze/bins_msg.csv", "w")
    for sample in ["Yap_AOA2", "Yap_AOA3", "Yap_AOA4"]:
        msg_title, MAGs_msg = collect_MAGs_msg(
            MAG_file_path.format(sample=sample),
            depth_file.format(sample=sample),
            checkM_file.format(sample=sample),
            iRep_file.format(sample=sample), gtdbtk_file.format(sample=sample))
