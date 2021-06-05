# -*- coding: utf-8 -*-
"""
 * @Editor: Lv Yongxin
 * @Date: 2020/08/29 11:08:50
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-06-05 21:49:16
 * @FilePath: /metaSC/PyLib/seqPipe/x03.2_integrate_KO.py
 * @Description: 对最多三种来源的KO注释结果整合

    usage: python x03.2_integrate_KO.py -p pattern [-o out_file]

    description: merge KO annotations from
        GhostKOALA, KofamKOALA, eggnog as the order
        Unfortunately, only those 3 can be recognized

    args:
        -p  pattern of the file
                'ghost', 'kofam' or 'eggnog' should in its filename
        -o  output file, default is stdout

"""


import os
import re
from sys import argv, stderr, stdout
from typing import Tuple

from PyLib.PyLibTool.file_info import verbose_import


verbose_import(__name__, __doc__)


def readline_GhostKOALA(line) -> Tuple[str, str]:
    """
        @description: read given file with Ghostkoala format
            and output gene -> ko map
    """
    temp = line.strip().split('\t')
    if len(temp) > 1:
        gene = temp[0]
        ko = temp[1]
        return gene, ko


def readline_KofamKOALA(line) -> Tuple[str, str]:
    """
        @description: read given file with KofamKOALA format
            and output gene -> ko map
    """
    temp = line.strip().split('\t')
    if len(temp) > 1:
        gene = temp[1]
        ko = temp[2]
        return gene, ko


def readline_eggnog(line) -> Tuple[str, str]:
    """
        @description: read given file with eggnog format
            and output gene -> ko map
    """
    if line[0] != '#':
        temp = line.strip().split('\t')
        ko = temp[8].split(',')[0][3:]
        if ko != '':
            gene = temp[0]
            return gene, ko


formats_map = {
    "ghost": readline_GhostKOALA,
    "kofam": readline_KofamKOALA,
    "eggnog": readline_eggnog,
}
formats_list = [
    "ghost", "kofam", "eggnog",
]


def main():
    # no input -> export help
    if len(argv) == 1:
        print(__doc__)
        print()
        exit(1)

    if len(argv) > 1:
        pattern = ".txt"
        in_dir = "./"
        out_name = ""
        i = 1
        while i < len(argv):

            if argv[i][0] == "-":  # key
                key = argv[i]
                if key == "-p":
                    pattern = argv[i + 1]
                    index = pattern.rfind("/") + 1
                    if index:
                        in_dir = pattern[:index]
                        pattern = pattern[index:]

                elif key == "-o":
                    out_name = argv[i + 1]
                i += 2

            else:  # what will happen then?
                i += 1
        stderr.write("pattern: {}\nin_dir: {}\nout_name: {}\n".format(
            pattern, in_dir, out_name))

        # find avaiable annotation file
        pattern = re.compile(pattern)
        ann_files = ["", "", ""]
        for file in sorted(os.listdir(in_dir)):
            if pattern.search(file):
                for i, source in enumerate(formats_list):
                    if source in file.lower():
                        ann_files[i] = os.path.join(in_dir, file)

        gene_KOs = {}
        for format, file in zip(formats_list, ann_files):
            if not file:
                stderr.write(format + " is not exist, skip\n")
                continue
            stderr.write("reading {}\n".format(file))
            with open(file) as fin:
                for line in fin:
                    hit = formats_map[format](line)
                    if hit and hit[0] not in gene_KOs:
                        gene_KOs[hit[0]] = hit[1]
            stderr.write("{} genes annotated...\n".format(len(gene_KOs)))

        if not gene_KOs:
            print("No KO found, pleace make sure there are"
                  " annotation files in your dictionary!")
            print(__doc__)
            exit(1)

        fout = open(out_name, "w") if out_name else stdout
        for gene, KO in gene_KOs.items():
            fout.write(gene + '\t' + KO + '\n')
        if out_name:
            fout.flush()
            fout.close()


if __name__ == "__main__":
    main()
    stderr.write("done!")
