# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-30 20:05:10
 * @Editors: Hwrn, LYX
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-12-13 17:12:54
 * @FilePath: /metaSC/PyLib/seqPipe/x03_link_gene.py
 * @Description:
    1.  it can generate/read subset gene or contig file, either fasta or table format.
    2.  it can intergrate KO
    3.  it can count gene and calculate reads hits/gene length (RPb)
        - Warning: RPb / mapped reads cannot reflact sample difference as well as tpm of gene only?
        - Modified: RPb_{gene} / sum(RPb_{markers}) can be used to compare between samples.
    4.  it can output a table of 'gene, KO, RPb, contig' (contig is for binning later)
"""

import argparse
import logging
import os
import re
from datetime import datetime
from io import FileIO
from sys import stdout
from typing import Dict, Hashable, Iterable, List, Set, Tuple

from Bio.SeqIO.FastaIO import SimpleFastaParser
from PyLib.reader.iters import emapper_iter, featureCounts_iter, read_table

logger = logging.getLogger(__name__)


def drop_gene(gene: str, subsets: Tuple[Hashable, bool]) -> bool:
    subset, subset_is_contig = subsets
    if subset:
        if subset_is_contig:
            return gene.rsplit('_', 1)[0] not in subset
        else:
            return gene not in subset
    return False


## collect gene KO
formats_list = [
    "ghost", "kofam", "eggnog",
]


def get_gene_KOs(ann_files: List[str],
                 subsets: Tuple[Hashable, bool] = None
                 ) -> Dict[str, str]:
    def GhostKOALA_iter(text: FileIO) -> Iterable[Tuple[str, str]]:
        for values in read_table(text):
            if len(values) > 1:
                gene, ko = values[0:2]
                yield gene, ko

    def KofamKOALA_iter(text: FileIO) -> Iterable[Tuple[str, str]]:
        for values in read_table(text):
            if len(values) > 1:
                gene, ko = values[1:3]
                yield gene, ko

    def eggnog_iter(text: FileIO) -> Iterable[Tuple[str, str]]:
        i_KEGG_ko = 11
        for values in emapper_iter(text):
            kos = values[i_KEGG_ko]
            if kos:
                gene = values[0]
                yield gene, kos.split(',')[0][3:]

    formats_func = {
        "ghost": GhostKOALA_iter,
        "kofam": KofamKOALA_iter,
        "eggnog": eggnog_iter,
    }

    gene_KOs: Dict[str, str] = {}
    for format, file in zip(formats_list, ann_files):
        if not os.path.exists(file):
            logger.warning(f"{format} file does not exist, skip")
            continue
        logger.info(f"reading {file}")

        with open(file) as file_in:
            for gene, ko in formats_func[format](file_in):
                if subsets and drop_gene(gene, subsets):
                    continue
                gene_KOs.setdefault(gene, ko)
        logger.warning(f"{len(gene_KOs)} genes annotated...")

    return gene_KOs


def get_gene_RPb(countfile, subsets
                 ) -> Tuple[List[str], Dict[str, List[float]]]:
    """
       @param {countlist: list}
            list read in by {callable: read_featureCounts}
       @return {samples: List, gene_RPb: dict} -> [samples, ], {
            contigName: {geneId: RPb}
        }
       @discussion: TPM is not avaiable only with featureCounts.
            In metagenome, totalMappedReads is at contig level
    """
    gene_RPb: Dict[str, Tuple[int, List[int]]] = {}

    logger.info(f"reading {countfile}")
    with open(countfile) as file_in:
        samples = next(featureCounts_iter(file_in))[6]

        for values in featureCounts_iter(file_in):
            geneId, Chr = values[:2]
            gene = f'{Chr}_{geneId.split("_")[1]}'

            if drop_gene(gene, subsets):
                continue
            length, readCounts = values[5:7]
            gene_RPb[gene] = [readCount / length for readCount in readCounts]

    logger.info(f"total {len(samples)} bam")
    return samples, gene_RPb


def main(subsets: Tuple[Hashable, bool],
         ann_files: List[str], countfile: str,
         show_contig: bool, output: FileIO):
    #logger.info(subsets)
    # 考虑到分支预测, 通过判断降低内存占用是更经济的
    gene_KOs = get_gene_KOs(
        ann_files, subsets) if any(ann_files) else {}
    samples, gene_RPb = get_gene_RPb(
        countfile, subsets) if countfile else ([], {})

    title, output_dict = ['#gene'], {gene: [] for gene
                                     in gene_KOs.keys() | gene_RPb.keys()}
    if gene_KOs:
        logger.info('collect KO messages')
        title += ['KO']
        for gene, output_list in output_dict.items():
            output_list.append(gene_KOs.get(gene, ''))
    if show_contig:
        logger.info('collect contig messages')
        title += ['contig']
        for gene, output_list in output_dict.items():
            output_list.append(gene.rsplit('_', 1)[0])
    if gene_RPb:
        logger.info('collect RPb messages')
        title += samples
        for gene, output_list in output_dict.items():
            output_list.extend(gene_RPb.get(gene, [''] * len(samples)))

    with output:
        print(*title, sep='\t', file=output)
        for gene in sorted(output_dict):
            print(gene, *output_dict[gene], sep='\t', file=output)
        logger.info(f'write to file {output.name}')

    return 0


def get_subset(raw_subset: List[str], threshold: int):
    subset = set()
    for file in raw_subset:
        total, keep = 0, 0
        suffix = file.rsplit('.', 1)[1]
        if suffix.startswith('f') and suffix.endswith('a'):
            # fasta format
            with open(file) as fasta_in:
                for title, seq in SimpleFastaParser(fasta_in):
                    total += 1
                    if threshold <= len(seq):
                        keep += 1
                        subset.add(title.split()[0])
        elif suffix.endswith('sv'):
            # table format
            sep = '\t' if suffix.startswith('t') else ','
            with open(file) as table_in:
                for line in read_table(table_in, sep=sep):
                    # the first column is title
                    total += 1
                    title = line[0]
                    keep += 1
                    subset.add(title)
        else:
            total += 1
            title = file
            keep += 1
            subset.add(title)
        logger.warning(f'{keep} of {total} featuers will be kept in {file}')
    return subset


def get_args() -> Tuple[Tuple[Set, bool],  # subset
                        List[str],  # KO
                        str,  # RPb
                        bool, FileIO]:
    parser = argparse.ArgumentParser(description=__doc__)
    set_args(parser)
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel.upper())  # info

    subsets = set(), False
    if args.subset:
        subset_threshold = args.subset_threshold
        if not isinstance(subset_threshold, int):
            logger.warning("subset_threshold must be an integer")
            subset_threshold = 0
            logger.warning(f"ignore parameter to {subset_threshold}")
        subset_feature_type = args.subset_feature_type
        if subset_feature_type == 'gene':
            subset_is_contig = False
        elif subset_feature_type == 'contig':
            subset_is_contig = True
        else:
            logger.warning('No other format, exit')
            parser.print_help()
            exit(1)

        subsets = get_subset(args.subset, subset_threshold), subset_is_contig

    ann_files = ["", "", ""]
    pattern = args.ko
    if pattern:
        in_dir = ""
        index = pattern.rfind("/") + 1
        if index:
            in_dir = pattern[:index]
            pattern = pattern[index:]
            if not os.path.exists(in_dir):
                logger.warning('illigal pattern, ignore')
                pattern = ''
        # find avaiable annotation file
        pattern = re.compile(pattern)
        for file in sorted(os.listdir(in_dir)):
            if pattern.search(file):
                for i, source in enumerate(formats_list):
                    if source in file.lower():
                        ann_files[i] = os.path.join(in_dir, file)

        if not any(ann_files):
            parser.print_help()
            logger.fatal('illigal pattern, please check!')
            exit(1)

    countfile = args.RPb
    if countfile:
        if not os.path.isfile(countfile):
            parser.print_help()
            logger.fatal('illigal countfile, please check!')
            exit(1)

    elif not any(countfile):
        parser.print_help()
        logger.fatal('either pattern or countfile should be provided!')
        exit(1)

    show_contig = args.show_contig

    output = args.output
    if output == 'stdout':
        output = stdout
    else:
        output = open(output, 'w')

    return (subsets,
            ann_files, countfile,
            show_contig, output)


def set_args(parser: argparse.ArgumentParser):
    parser.add_argument('--loglevel', default='INFO', type=str,
                        help='set level of logger')

    parser.add_argument('-s', '--subset', nargs='*', type=str,
                        help='Subset files to reduce geneset size. '
                             'If subset is specified, then only genes match '
                             '"--subset" will be recorded. '
                             'Automatically recorded the suffix of files: '
                             'should be in ["f*a", "tsv", "csv"]; '
                             'or just regarded as feature names'
                        )
    parser.add_argument('--subset-threshold', default=33, type=int,
                        help='Threshold for feature to keep if supported. '
                             'Only useful when "--subset" is specified. '
                             'Only feature is shorter then this threshold '
                             'will be discarded. '
                             '\n    e.g. If 33 <= len(gene A), '
                             'then gene A will be kept. '
                             'Else if set to 0 or no matched features in '
                             '"--subset", all genes will be kept.'
                        )
    parser.add_argument('--subset-feature-type', default='gene',
                        help='The feature type of subset. '
                             'Only useful when "--subset" is specified. '
                             'Only ["gene", "contig"] is supported. '
                             'If "--subset" is contig, prefix of gene will '
                             'be checked.'
                        )

    parser.add_argument('--ko', '--KO-pattern', default='', type=str,
                        help='Pattern of KO files. KO annotations from '
                             '"GhostKOALA, KofamKOALA, eggnog" '
                             '(just in this order) '
                             'will be merged and output. '
                             '"ghost", "kofam" or "eggnog" MUST in its filename.'
                        )
    parser.add_argument('--RPb', '--featuer-counts', default='', type=str,
                        help='File generated by featureCounts '
                             'for calulating TPM. '
                             '> TPM is used to show the rate of a feature to be '
                             'recognized at this depth '
                             'if gene length is the same. '
                             'TPM = transcript / total_transcript * 1000000. '
                             'e.g. "featureCounts '
                             '    -a ${gff} '
                             '    -o <in_file> '
                             '    -t CDS -g ID '
                             '    -p ${bam}"'
                        )
    parser.add_argument('--show-contig', action='store_true',
                        help='Keep a new line of contig name. '
                             'It can be useful if someone want to select gene '
                             'according to contigs.'
                        )
    parser.add_argument('-o', '--output', type=str, default='stdout',
                        help='output file.')


def run():
    args = get_args()

    now = datetime.now()
    logger.warning('>>> job start at ' + now.strftime('%Y-%m-%d %H:%M:%S'))
    state = main(*args)
    logger.warning('>>> job run time: ' + str(datetime.now() - now))
    if state == 0:
        logger.info('success!')


if __name__ == '__main__':
    run()
