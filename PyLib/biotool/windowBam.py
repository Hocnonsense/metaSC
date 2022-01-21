# -*- coding: utf-8 -*-
"""
 * @Date: 2022-01-21 15:59:44
 * @Editors: derek.bickhart-adm, Hwrn
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-21 16:27:33
 * @FilePath: /metaSC/PyLib/biotool/windowBam.py
 * @Description:
"""

from typing import Dict, List
import pysam
import click
from collections import defaultdict


def windowBam(bam, output, windowlength):
    with pysam.AlignmentFile(bam, 'rb') as bamfile:
        # create windows
        references = bamfile.references
        lengths = bamfile.lengths
        winlist: Dict[str, List[window]] = defaultdict(list)
        for r, l in zip(references, lengths):
            for i in range(0, l, windowlength):
                winlist[r].append(window(r, i, i + windowlength))

        for c, w in winlist.items():
            for i, win in enumerate(w):
                count = 0
                mapq0 = 0
                for s in bamfile.fetch(c, win.start, win.end):
                    if s.is_secondary:
                        continue
                    count += 1
                    if s.mapping_quality == 0:
                        mapq0 += 1
                winlist[c][i].recordCount(mapq0, count)

    # Now print it all out
    with open(output, 'w') as out:
        for c, w in winlist.items():
            for win in w:
                print(win.getBed(), file=out)

class window:

    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start
        self.end = end
        self.mapq0 = 0
        self.totReadCount = 0

    def recordCount(self, mapq0, totReads):
        self.mapq0 = mapq0
        self.totReadCount = totReads

    def getRatio(self):
        return self.mapq0 / self.totReadCount

    def hasMapq0(self):
        return self.mapq0 != 0

    def getBed(self):
        return "\t".join([self.contig, self.start, self.end, self.hasMapq0, self.getRatio, self.totReadCount, self.mapq0])


@click.command()
@click.option("-b", "--bam", required=True, type=str,
              type=click.Path(exists=True, dir_okay=False),
              help="Input bam file. Must be coord sorted and indexed.")
@click.option('-o', '--output', required=True, type=str,
              help="Output file name. Output format is bed format.")
@click.option('-w', '--windowlength', default = 5000, type=int,
              help="Window length for contig windows")
def arg_parse(bam, output, windowlength):
    """Window analysis to calculate ratio of mapq0 reads in bam files"""
    return bam, output, windowlength


def run():
    windowBam(*arg_parse())


if __name__ == "__main__":
    run()
