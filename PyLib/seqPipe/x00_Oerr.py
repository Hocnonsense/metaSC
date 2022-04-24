# -*- coding: utf-8 -*-
"""
 * @Date: 2021-07-01 20:30:00
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-24 19:03:23
 * @FilePath: /metaSC/PyLib/seqPipe/x00_Oerr.py
 * @Description:
"""

import os
import re
from ast import literal_eval as eval
import sys
from typing import Callable, Dict, List, Union, TextIO

from PyLib.reader.read_outputs import fasta, jgi_depths
from PyLib.biotool.fna_msg import statistic_fna, seq_total_depth


line = ""


def sickle_mapper(text: TextIO):
    """{
        'header': ['key',
                   'input reads',
                   'paired reads kept', 'SE1 reads kept', 'SE2 reads kept',
                   'paired reads discarded', 'SE1 reads discarded', 'SE2 reads discarded'],
        'demo': '''
            #2021-05-28 20:58:21

            PE forward file: /lustre/home/acct-clsxx/clsxx/Data/Database2/BaiduNetdiskDownload/YW.020/YW.020_1.fq.gz
            PE reverse file: /lustre/home/acct-clsxx/clsxx/Data/Database2/BaiduNetdiskDownload/YW.020/YW.020_2.fq.gz

            Total input FastQ records: 118472086 (59236043 pairs)

            FastQ paired records kept: 117609902 (58804951 pairs)
            FastQ single records kept: 386404 (from PE1: 238384, from PE2: 148020)
            FastQ paired records discarded: 89376 (44688 pairs)
            FastQ single records discarded: 386404 (from PE1: 148020, from PE2: 238384)
            ''',
        'suffix': '.out'
    }"""
    global line
    for line in text:
        if line.startswith("PE forward file:"):
            break
    key = line.strip()

    for line in text:
        if line.startswith("Total input FastQ records:"):
            break
    in_bases, _ = re.match(  # type: ignore
        re.compile(r"^Total input FastQ records: (\d+) \((\d+) pairs\)$"), line
    ).groups()
    line = text.readline()

    line = text.readline()
    pair_keep, _ = re.match(  # type: ignore
        re.compile(r"^FastQ paired records kept: (\d+) \((\d+) pairs\)$"), line
    ).groups()
    line = text.readline()
    _, s1_keep, s2_keep = re.match(  # type: ignore
        re.compile(
            r"^FastQ single records kept: (\d+) \(from PE1: (\d+), from PE2: (\d+)\)$"
        ),
        line,
    ).groups()
    line = text.readline()
    pair_discard, _ = re.match(  # type: ignore
        re.compile(r"^FastQ paired records discarded: (\d+) \((\d+) pairs\)$"), line
    ).groups()
    line = text.readline()
    _, s1_discard, s2_discard = re.match(  # type: ignore
        re.compile(
            r"^FastQ single records discarded: (\d+) \(from PE1: (\d+), from PE2: (\d+)\)$"
        ),
        line,
    ).groups()

    return (
        key,
        in_bases,
        pair_keep,
        s1_keep,
        s2_keep,
        pair_discard,
        s1_discard,
        s2_discard,
    )


def bbmap_mapper(text: TextIO, lastline: str = None):
    """{
        'header': scaffold, [
            'Mapped Read 1 File',
            'Reads', 'Mapped reads', 'Mapped bases',
            'Percent mapped', 'Percent proper pairs',
            'Average coverage', 'Average coverage with deletions',
            'Standard deviation', 'Percent of reference bases covered'
        ]
    }"""
    if lastline is None:
        for lastline in text:
            match = re.match(  # type: ignore
                re.compile(
                    r"^Executing align2.BBMap \[.+, in=([^\s]+), .*ref=([^\s]+), "
                ),
                lastline,
            )
            if match:
                break
    match = re.match(  # type: ignore
        re.compile(r"^Executing align2.BBMap \[.+, in=([^\s]+), .*ref=([^\s]+), "),
        lastline,
    )
    if not match:
        return
    (in1_trim, scf) = match.groups()

    for match in (
        re.match(  # type: ignore
            re.compile(r"^Reads:                               	(\d+)$"), line
        )
        for line in text
    ):
        if match:
            (Reads,) = match.groups()
            break
    line = text.readline()
    (M_reads,) = re.match(  # type: ignore
        re.compile(r"^Mapped reads:                        	(\d+)$"), line
    ).groups()
    line = text.readline()
    (M_bases,) = re.match(  # type: ignore
        re.compile(r"^Mapped bases:                        	(\d+)$"), line
    ).groups()
    line = text.readline()
    # re.match(re.compile(r'^Ref scaffolds:                       	(\d+)$'), line).groups()[0]
    line = text.readline()
    # re.match(re.compile(r'^Ref bases:                           	(\d+)$'), line).groups()[0]
    line = text.readline()

    line = text.readline()
    (P_mapped,) = re.match(  # type: ignore
        re.compile(r"^Percent mapped:                      	(\d+.\d+)$"), line
    ).groups()
    line = text.readline()
    (P_pairs,) = re.match(  # type: ignore
        re.compile(r"^Percent proper pairs:                	(\d+.\d+)$"), line
    ).groups()
    line = text.readline()
    (Avg_cover,) = re.match(  # type: ignore
        re.compile(r"^Average coverage:                    	(\d+.\d+)$"), line
    ).groups()
    line = text.readline()
    (Avg_cov_del,) = re.match(  # type: ignore
        re.compile(r"^Average coverage with deletions:     	(\d+.\d+)$"), line
    ).groups()
    line = text.readline()
    (SD,) = re.match(  # type: ignore
        re.compile(r"^Standard deviation:                    	(\d+.\d+)$"), line
    ).groups()
    line = text.readline()
    # re.match(re.compile(r'^Percent scaffolds with any coverage: 	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    (P_covered,) = re.match(  # type: ignore
        re.compile(r"^Percent of reference bases covered:  	(\d+.\d+)$"), line
    ).groups()
    return scf, (
        in1_trim,
        Reads,
        M_reads,
        M_bases,
        P_mapped,
        P_pairs,
        Avg_cover,
        Avg_cov_del,
        SD,
        P_covered,
    )


def assem_mapper(text: TextIO):
    """{
        'header': [
            'key',
            'SeqNumbers', 'MaxLength', 'GenomeSize',
            'GC', 'N50', 'L50',
            'total_depth',
            'Mapped Read 1 File,',
            'Reads,', 'Mapped reads,', 'Mapped bases,',
            'Percent mapped,', 'Percent proper pairs,',
            'Average coverage,', 'Average coverage with deletions,',
            'Standard deviation,', 'Percent of reference bases covered,'
        ],
        'demo': '''
            Reads:                               	124372186
            Mapped reads:                        	116818023
            Mapped bases:                        	17630289518
            Ref scaffolds:                       	223406
            Ref bases:                           	375304290

            Percent mapped:                      	93.926
            Percent proper pairs:                	90.296
            Average coverage:                    	46.976
            Average coverage with deletions:     	46.827
            Standard deviation:                    	179.611
            Percent scaffolds with any coverage: 	100.00
            Percent of reference bases covered:  	99.75
            ''',
        'suffix': '.err',
        'desc': '''
            some times, multipile bbmap will be called, and should handle this
        '''
    }"""
    global line
    bbmap_outs = []
    for line in text:
        bbmap_out = bbmap_mapper(text, line)
        if bbmap_out:
            bbmap_outs.append(bbmap_out[1])
            scf = bbmap_out[0]  # only keep the last result
        match = re.match(re.compile(r"^Output depth matrix to ([^\s]+)$"), line)
        if match:
            (depth,) = match.groups()
            break

    try:
        fna_msg = statistic_fna(fasta(scf))
        with open(depth) as fi:
            sample_list, ctg_depth = jgi_depths(fi)
        totalAvgDepth, depths = seq_total_depth(ctg_depth, fasta(scf))
        values = list(fna_msg) + [totalAvgDepth]
    except FileNotFoundError:
        values = [""] * (6 + 1)

    return (
        scf,
        *values,
        *[",".join(i) for i in zip(*bbmap_outs)],
    )


def collect_folder(path: str, mapper: Callable[[TextIO], List[str]], suffix=""):
    results = {}
    if not suffix:
        suffix = eval(mapper.__doc__)["suffix"]  # type: ignore
    files = os.listdir(path)
    for file in (os.path.join(path, _file) for _file in sorted(files)):
        if file.endswith(suffix):
            with open(file) as file_in:
                try:
                    values = mapper(file_in)
                    results[values[0]] = values[1:]
                except AttributeError as e:
                    print(f"AttributeError at file {file}", file=sys.stderr)
                    print("error line:\n", line, file=sys.stderr)
                except Exception as e:
                    print(e, file=sys.stderr)
                    print("at file:", file, file=sys.stderr)
                    print("line:\n", line, file=sys.stderr)
    return eval(mapper.__doc__)["header"], results  # type: ignore


def print_table(
    header: Union[List, Dict],
    results: Union[List, Dict] = None,
    output: TextIO = sys.stdout,
    sep="\t",
):
    if not results:
        results = header
        header = []
    if header:
        print("#", end="", file=output)
        print(*header, sep=sep, file=output)
    if isinstance(results, dict):
        for key, values in results.items():
            print(key, *values, sep=sep, file=output)
    elif isinstance(results, list):
        for values in results:
            print(values, sep=sep, file=output)


def test():
    print_table(*collect_folder("log/01.2_trim", sickle_mapper), output=sys.stdout)
    print_table(*collect_folder("log/02.3_depth", assem_mapper), output=sys.stdout)


if __name__ == "__main__":
    test()
