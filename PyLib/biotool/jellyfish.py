# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-16 11:13:39
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-07-30 19:41:02
 * @filePath: /2021_Spr-Cources/Omics/final/p1/02.py
 * @Description:
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def KMER_TABLE_FORMAT(kmer):
    return f'{kmer}_count.tsv'


def set_backgroud_frequency(filepath: str):
    # assume file exists
    base_backgroud_count = {base: 0.0 for base in 'ATGC'}
    with open(filepath) as fa_in:
        for line in fa_in:
            if line.startswith('>'):
                logger.info(line)
            else:
                line = line.upper()
                for base in base_backgroud_count:
                    base_backgroud_count[base] += line.count(base)

    total = sum(base_backgroud_count.values())
    logger.info(base_backgroud_count)

    return {base: count / total for (base, count)
            in base_backgroud_count.items()}


def expected_frequency(sequence: str, base_backgroud_frequency: Dict[str, float]):
    sequence = sequence.upper()
    frequency = 1.0
    for base, base_frequency in base_backgroud_frequency.items():
        frequency *= base_frequency ** sequence.count(base)
    return frequency


def most_n_kmers(kmer, kmer_len, base_backgroud_frequency):
    def min_kmer(kmers: Tuple[List[str], List[int]]):
        i = min(range(kmer_len), key=lambda x: kmers[1][x])
        return i, kmers[1][i]
    most_freq_kmers = [''] * kmer_len, [0.0] * kmer_len
    most_rich_kmers = [''] * kmer_len, [0.0] * kmer_len
    ## run jellyfish
    """
    jellyfish count \
        -s $MEMORY -t $THREAD -C 00_data/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna \
        -m $k -o "$k"_count.jf

    jellyfish dump \
        -c -t \
        "$k"_count.jf \
    > {KMER_TABLE_FORMAT('$k')}
    """
    with open(KMER_TABLE_FORMAT(kmer)) as table_in:
        table_iter = ((kseq, int(count)) for (kseq, count) in
                      (line.split() for line in table_in))
        for i, (kseq, count) in zip(range(kmer_len), table_iter):
            most_freq_kmers[0][i] = kseq
            most_freq_kmers[1][i] = count
            most_rich_kmers[0][i] = kseq
            most_rich_kmers[1][i] = count / expected_frequency(
                kseq, base_backgroud_frequency)
        p_min_freq, min_freq = min_kmer(most_freq_kmers)
        p_min_rich, min_rich = min_kmer(most_rich_kmers)
        for kseq, count in table_iter:
            if min_freq < count:
                most_freq_kmers[0][p_min_freq] = kseq
                most_freq_kmers[1][p_min_freq] = count
                # update
                p_min_freq, min_freq = min_kmer(most_freq_kmers)
            if min_rich < count:
                most_rich_kmers[0][p_min_rich] = kseq
                most_rich_kmers[1][p_min_rich] = count
                # update
                p_min_rich, min_rich = min_kmer(most_rich_kmers)
    return ([most_freq_kmers[0][i]
             for i in sorted(range(kmer_len),
                             key=lambda i: most_freq_kmers[1][i], reverse=True)],
            [most_rich_kmers[0][i]
             for i in sorted(range(kmer_len),
                             key=lambda i: most_rich_kmers[1][i], reverse=True)])


def main():
    base_backgroud_frequency = set_backgroud_frequency(
        '00_data/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna')
    with open('most_freq_kmers.tsv', 'w') as freq_out, \
            open('most_rich_kmers.tsv', 'w') as rich_out:
        print('kmer', range(3, 21), sep='\t', file=freq_out)
        print('kmer', range(3, 21), sep='\t', file=rich_out)
        for kmer in range(3, 21):
            most_freq_kmers, most_rich_kmers = most_n_kmers(
                kmer, 20, base_backgroud_frequency)
            print(kmer, *most_freq_kmers, sep='\t', file=freq_out)
            print(kmer, *most_rich_kmers, sep='\t', file=rich_out)


if __name__ == '__main__':
    main()
