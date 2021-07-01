# -*- coding: utf-8 -*-
"""
 * @Date: 2021-07-01 20:30:00
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-07-01 23:31:42
 * @FilePath: /metaSC/PyLib/seqPipe/x00_Oerr.py
 * @Description:
"""

import os
import re
from ast import literal_eval as eval
from io import FileIO
import sys
from typing import Callable, Dict, List, Union

from PyLib.reader.read_outputs import fasta, jgi_depths
from PyLib.biotool.fna_msg import statistic_fna, seq_total_depth


def sickle_mapper(text: FileIO):
    """{
        'header': ['key', 'input reads',
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
    for line in text:
        if line.startswith('PE forward file:'):
            break
    key = line.strip()

    for line in text:
        if line.startswith('Total input FastQ records:'):
            break
    in_bases, _ = re.match(re.compile(r'^Total input FastQ records: (\d+) \((\d+) pairs\)$'), line).groups()
    line = text.readline()

    line = text.readline()
    pair_keep, _ = re.match(re.compile(r'^FastQ paired records kept: (\d+) \((\d+) pairs\)$'), line).groups()
    line = text.readline()
    _, s1_keep, s2_keep = re.match(re.compile(r'^FastQ single records kept: (\d+) \(from PE1: (\d+), from PE2: (\d+)\)$'), line).groups()
    line = text.readline()
    pair_discard, _ = re.match(re.compile(r'^FastQ paired records discarded: (\d+) \((\d+) pairs\)$'), line).groups()
    line = text.readline()
    _, s1_discard, s2_discard = re.match(re.compile(r'^FastQ single records discarded: (\d+) \(from PE1: (\d+), from PE2: (\d+)\)$'), line).groups()

    return key, in_bases, pair_keep, s1_keep, s2_keep, pair_discard, s1_discard, s2_discard


def assem_mapper(text: FileIO):
    """{
        'header': [
            'SeqNumbers', 'MaxLength', 'GenomeSize',
            'GC', 'N50', 'L50',
            'total_depth', '{depths}',
            'Reads', 'Mapped reads', 'Mapped bases',
            'Percent mapped', 'Percent proper pairs',
            'Average coverage', 'Average coverage with deletions', 'Standard deviation'
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
        'suffix': '.err'
    }"""
    for line in text:
        if line.startswith('+ scf='):
            break
    scf = re.match(re.compile(r'\+ scf=(.+)$'), line).groups()[0]

    for line in text:
        if line.startswith('+ depth='):
            break
    depth = re.match(re.compile(r'\+ depth=(.+)$'), line).groups()[0]
    try:
        fna_msg = statistic_fna(fasta(scf))
        with open(depth) as fi:
            sample_list, ctg_depth = jgi_depths(fi)
        totalAvgDepth, depths = seq_total_depth(ctg_depth, fasta(scf))
        values = list(fna_msg) + [totalAvgDepth] + depths
    except FileNotFoundError:
        values = []

    for line in text:
        if line.startswith('   ------------------   Results   ------------------'):
            break
    for line in text:
        if line.startswith('Reads:                               	'):
            break
    Reads = re.match(re.compile(r'^Reads:                               	(\d+)$'), line).groups()[0]
    line = text.readline()
    M_reads = re.match(re.compile(r'^Mapped reads:                        	(\d+)$'), line).groups()[0]
    line = text.readline()
    M_bases = re.match(re.compile(r'^Mapped bases:                        	(\d+)$'), line).groups()[0]
    line = text.readline()
    # re.match(re.compile(r'^Ref scaffolds:                       	(\d+)$'), line).groups()[0]
    line = text.readline()
    # re.match(re.compile(r'^Ref bases:                           	(\d+)$'), line).groups()[0]
    line = text.readline()

    line = text.readline()
    P_mapped = re.match(re.compile(r'^Percent mapped:                      	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    P_pairs = re.match(re.compile(r'^Percent proper pairs:                	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    Avg_cover = re.match(re.compile(r'^Average coverage:                    	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    Avg_cov_del = re.match(re.compile(r'^Average coverage with deletions:     	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    SD = re.match(re.compile(r'^Standard deviation:                    	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    #re.match(re.compile(r'^Percent scaffolds with any coverage: 	(\d+.\d+)$'), line).groups()[0]
    line = text.readline()
    P_covered = re.match(re.compile(r'^Percent of reference bases covered:  	(\d+.\d+)$'), line).groups()[0]

    return [scf] + values + [Reads, M_reads, M_bases,
                             P_mapped, P_pairs,
                             Avg_cover, Avg_cov_del,
                             SD, P_covered]


def collect_folder(path: str, mapper: Callable[[FileIO], List[str]], suffix=''):
    results = {}
    if not suffix:
        suffix = eval(mapper.__doc__)['suffix']
    files = os.listdir(path)
    for file in (os.path.join(path, _file) for _file in sorted(files)):
        if file.endswith(suffix):
            with open(file) as file_in:
                try:
                    values = mapper(file_in)
                    results[values[0]] = values[1:]
                except AttributeError as e:
                    print(f'AttributeError at file {file}')
                    print(e)
    return eval(mapper.__doc__)['header'], results


def print_table(output: FileIO,
                header: List = [], results: Union[List, Dict] = [],
                sep='\t'):
    if not results:
        results = header
        header = []
    if header:
        print('#', end='', file=output)
        print(*header, sep=sep, file=output)
    if isinstance(results, dict):
        for key, values in results.items():
            print(key, *values, sep=sep, file=output)
    elif isinstance(results, list):
        for values in results:
            print(values, sep=sep, file=output)


def test():
    print_table(sys.stdout, *collect_folder('Oerr/01.2_trim', sickle_mapper))
    print_table(sys.stdout, *collect_folder('Oerr/02.3_depth', assem_mapper))


if __name__ == '__main__':
    test()
