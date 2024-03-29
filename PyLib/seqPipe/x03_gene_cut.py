# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-24 10:24:10
 * @Editor: LYX
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-08-08 20:59:53
 * @FilePath: /metaSC/PyLib/seqPipe/x03_gene_cut.py
 * @Description:
        update from LYX's script
    x03.1_gene_cut.py <in_file_prefix> <out_file_prefix> <threshold> [-h] [--help]
        @description:   Remove sequence shorter than given threshold.
                        This file is designed to control quality of
                        prodigal.faa, prodigal.fna, prodigal.gff generated by
                        prodigal.
        @param:         <in_file_prefix>:   input file prefix. For example,
                                        if in_file_prefix is "./M001-prodigal",
                                        then will read these three files:
                                                ./M001-prodigal.faa
                                                ./M001-prodigal.fna
                                                ./M001-prodigal.gff
                                        **NOTE**: '"in_file_prefix".faa' must
                                        exist.
                        * <in_file_prefix>: other input file prefix.
                        <out_file_prefix>:  output file prefix. Similar to
                                        <in_file_prefix>
                        <threshold>:        threshold for aa, typtically, 33
                                        Any sequence in in_file_prefix.faa but
                                        < threshold will be discarded.
                                        Its record in in_file_prefix.fna and
                                        in_file_prefix.gff will be discarded as
                                        well.
"""

import os
from sys import argv, stderr
from Bio import SeqIO, SeqRecord
from typing import Iterable, List, Optional, Set

from PyLib.PyLibTool.file_info import verbose_import


verbose_import(__name__, __doc__)


class PrefixFiles:
    def __init__(self, prefix) -> None:
        self.prefix = os.path.abspath(os.path.expanduser(prefix))

    def __str__(self):
        return f"<prefix_files: {self.prefix}>"

    def __call__(self, suffix: str, check_exist: bool = None) -> str:
        """
        @param {check_exist: bool}:
             True: return False if file not exist
             False: return False if file exist
             None: do not check
        """
        file_name = self.prefix + suffix

        if os.path.isdir(file_name):
            raise IsADirectoryError
        if check_exist is None or check_exist is os.path.isfile(file_name):
            return file_name

        return ""


def fa_filt_iter(
    fas: Iterable[SeqRecord.SeqRecord],
    keepIds: Set[str],
    threshold: Optional[int] = None,
):
    """Any sequence will be kept if
    1.  in keepIds
    2.  not less than threshold

    @paran threshold: typtically, 33
    """
    keep_lens: List[int] = []
    discard_lens: List[int] = []
    for record in fas:
        if threshold:
            if len(record.seq) >= threshold:
                keepIds.add(record.id)
                keep_lens.append(len(record.seq))
            else:
                discard_lens.append(len(record.seq))
        if record.id in keepIds:
            yield record
    if discard_lens:
        print(f"discarded: {len(discard_lens)}, {sum(discard_lens)} bases(aa)")
        print(f"kept: {len(keep_lens)}, {sum(keep_lens)} bases(aa)")


def gff_filt_iter(gffin: Iterable[str], keepIds: Set[str]):
    for line in gffin:
        if line[0] == "#":
            continue
        temp = line.strip().split("\t")
        scaffold = temp[0]
        s_count = temp[8].split(";")[0].split("_")[1]
        sid = scaffold + "_" + s_count
        if sid in keepIds:
            yield line


def parse_args():
    if "-h" in argv or "--help" in argv or len(argv) == 1:
        exit(0)

    sc = argv[0]
    in_file_prefixs = [PrefixFiles(in_file_prefix) for in_file_prefix in argv[1:-2]]
    out_file_prefix = PrefixFiles(argv[-2])
    num = int(argv[-1])
    args = in_file_prefixs, out_file_prefix, num
    print(sc, *args, sep="\n" + " " * 4, file=stderr)
    return args


def main(in_file_prefix, num, out):
    keepIds = set()
    faao, fnao, gffo = out

    suffix = ".faa"
    in_file = in_file_prefix(suffix, True)
    print("read", in_file, file=stderr)
    if in_file:
        with open(in_file) as fin:
            for record in fa_filt_iter(SeqIO.parse(fin, "fasta"), keepIds, num):
                print(record.format("fasta"), file=faao)
    else:
        print("invalid", suffix, "file:", in_file, file=stderr)
        return

    suffix = ".fna"
    in_file = in_file_prefix(suffix, True)
    print("read", in_file, file=stderr)
    if in_file:
        with open(in_file) as fin:
            for record in fa_filt_iter(SeqIO.parse(fin, "fasta"), keepIds):
                print(record.format("fasta"), file=fnao)
    else:
        print(f"    {in_file_prefix(suffix)} no found, pass")

    suffix = ".gff"
    in_file = in_file_prefix(suffix, True)
    print("read", in_file, file=stderr)
    if in_file:
        with open(in_file) as fin:
            for line in gff_filt_iter(fin, keepIds):
                gffo.write(line)
    else:
        print("    {in_file} no found, pass".format(in_file=in_file_prefix(suffix)))


if __name__ == "__main__":
    in_file_prefixs, out_file_prefix, num = parse_args()

    with open(out_file_prefix(".faa", "w")) as faao, open(
        out_file_prefix(".fna", "w")
    ) as fnao, open(out_file_prefix(".gff", "w")) as gffo:
        for in_file_prefix in in_file_prefixs:
            main(in_file_prefix, num, (faao, fnao, gffo))
