# -*- coding: utf-8 -*-
"""
 * @Date: 2021-08-09 17:43:13
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-08-20 16:26:46
 * @FilePath: /metaSC/PyLib/biotool/featureReads.py
 * @Description:
"""

import os
import sys
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Tuple

from Bio import SeqFeature, SeqIO
from PyLib.PyLibTool.file_info import verbose_import

logger = verbose_import(__name__, __doc__)


try:
    import pysam
    from BCBio import GFF
except ImportError:
    logger.error('module pysam or BCBio is not installed, skipping')


def extract_exons(seq: SeqFeature.SeqFeature) -> List[SeqFeature.SeqFeature]:
    """
     * @description:
     * @param {SeqRecord} seq
     * @return {List[SeqFeature.SeqFeature]} List of SeqFeature at exon level
        region of Seqfeature: 0 ... [start, end)
        region of gff       : 1 ... [start, end]
        e.g.: exon:ENST00000620521.1:1
            Seqfeature: 64327417, 64327972, -1
            gff       : 64327418, 64327972, -

    """
    features: List[SeqFeature.SeqFeature] = seq.features.copy()
    exons: List[SeqFeature.SeqFeature] = []

    while features:
        feature = features.pop()
        if feature.type != 'exon':
            features.extend(feature.sub_features)
        else:
            start, end = feature.location.start.position, feature.location.end.position
            assert start < end
            exons.append(feature)
            #exon_seq :SeqFeature.SeqFeature = feature.extract(chr20)
            #exon_seq :Seq.Seq = exon_seq.seq if feature.strand == -1 else exon_seq.seq.reverse_complement()
            #exons.append(SeqFeature.SeqFeature(exon_seq, feature.id, feature.id, features=[feature]))
    return exons


def line_exons(exons: List[SeqFeature.SeqFeature]):
    """
     * @description:
        make a list of number that will never overlap
        it will start with 0 and any loc between [2k+1, 2k+2] is in a exon
        if multi exons at the same loc, goto exon_loc[2k] for more infomation
     * @param {*}
     * @return {*}
    """
    exon_sort: Callable[[SeqFeature.SeqFeature], Tuple[int, int]] = lambda exon: (exon.location.start.position, exon.location.end.position)

    exons.sort(key=exon_sort)
    locats: List[int] = [0]
    exon_loc: Dict[int, List[SeqFeature.SeqFeature]] = {}
    for exon in exons:
        start, end = exon.location.start.position, exon.location.end.position

        if start < locats[-1]:
            locats[-1] = end
        else:
            locats.extend([exon.location.start.position, exon.location.end.position])
        exon_loc.setdefault(len(locats) - 2, []).append(exon)
    #duplicates = {1: 2988, 2: 1392, 3: 957, 7: 241, 6: 356, 8: 181, 9: 184, 10: 156, 4: 586, 12: 55, 13: 53, 5: 485, 27: 2, 91: 1, 30: 5, 31: 3, 34: 6, 36: 2, 32: 2, 38: 3, 35: 3, 11: 89, 20: 21, 28: 2, 17: 21, 15: 39, 73: 1, 89: 1, 55: 2, 16: 32, 14: 24, 18: 32, 19: 24, 22: 10, 29: 6, 41: 1, 46: 2, 50: 3, 131: 1, 128: 1, 124: 1, 106: 1, 75: 1, 23: 6, 24: 7, 33: 2, 21: 9, 26: 4, 88: 1, 86: 1, 54: 1, 115: 1, 52: 2, 40: 1, 44: 2, 59: 1, 64: 1, 45: 2, 69: 1, 25: 5, 39: 3, 53: 2, 57: 3, 49: 1, 43: 3, 42: 1, 90: 1}
    #for _exons in exon_loc.values():
    #    duplicates[len(_exons)] = duplicates.get(len(_exons), 0) + 1
    return locats, exon_loc


def find_exons_left(locat: int, locats: List[int], index: Optional[Tuple[int, int]] = None) -> int:
    """
     * @description: find exons from seqs
     * @param {int}
        locat the location of query base
        must > 0
     * @param {List[int]} locats
        generated by line_exons
     * @param {Dict[int, ]} exon_loc
        generated by line_exons
     * @param {Optional} index
        [left, right) to cut the loacts
        geneted when running, please not change it
     * @return {*}
    """
    # from exon
    try:
        l, r = index
    except TypeError:
        l, r = 0, len(locats) - 1
        if locats[r] < locat:
            # ! bigger than max
            return
            # else all befpre
    # ! exist at left value
    if l + 1 >= r:
        if l % 2:  # 1, is the start of gene
            return l
        return
    mid = (l + r) // 2
    if locat < locats[mid]:
        return find_exons_left(locat, locats, (l, mid))
    #elif locat == locats[mid]:
    #    return mid
    else:
        return find_exons_left(locat, locats, (mid, r))


def exon_iter(locat: int, locats: List[int], exon_loc: Dict[int, List[SeqFeature.SeqFeature]]) -> Iterable[SeqFeature.SeqFeature]:
    left = find_exons_left(locat, locats)
    for exon in exon_loc.get(left, []):
        # ! [left, right)
        if exon.location.start.position <= locat and locat < exon.location.end.position:
            yield exon


def filter_OmicsPrj1(read: pysam.AlignedSegment) -> Iterator[pysam.AlignedSegment]:
    if read.mapping_quality < 2:
        yield read
    #mismatch = 0
    #for signal, seqlen in read.cigartuples:
    #    if signal not in {0,3}:
    #        mismatch += seqlen
    #if 2 < mismatch:
    #    return True


def filter_cmseq(read: pysam.AlignedSegment,
                 minlen=70, minqual=30, maxsnps=1.0, exclude_seqs: Dict[str]=None
                 ) -> Iterator[pysam.AlignedSegment]:
	alignment_len = int(read.query_alignment_length)

	qualities = read.query_qualities
	snps_rate =float(read.get_tag('NM')) / alignment_len
	refname = read.reference_name
	#meanqualities =np.mean(read.query_qualities)

	if (not read.is_secondary) \
            and (alignment_len >= minlen) \
            and (qualities >= minqual) \
            and (snps_rate <= maxsnps) \
            and not (exclude_seqs and refname in exclude_seqs):
		yield read


def load_bam_2_memory(filename: str) -> Dict[str, pysam.AlignedSegment]:
    uniq_bamis: Dict[str, List[pysam.AlignedSegment]] = {}
    with pysam.AlignmentFile(filename, 'rb') as bami:  # pysam.AlignmentFile
        for read in filter_OmicsPrj1(bami.fetch()):
            uniq_bamis.setdefault(read.query_name + ('_1' if read.is_read1 else '_2'), []).append(read)

    # Believe me duplicate reads don't exist
    assert not {tuple(reads) for reads in uniq_bamis.values() if len(reads) > 1}, \
           'some reads are mapped to multipile locats. please check'
    uniq_bamis = {read_name: read[0] for read_name, read in uniq_bamis.items()}

    return uniq_bamis


def read_junc(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """
    get --]...[--
    """
    lend, rstart = read.reference_start, read.reference_end
    beforeN = True
    for signal, seqlen in read.cigartuples:
        if signal == 3:
            beforeN = False
        elif signal == 0:
            if beforeN:
                lend += seqlen
            else:
                rstart -= seqlen
    return lend, rstart


def map_read_exon(uniq_bamis: Dict[str, pysam.AlignedSegment],
                  locats: List[int],
                  exon_loc: Dict[int, List[SeqFeature.SeqFeature]],
                  anno_method: str
                  ) -> Tuple[Dict[SeqFeature.SeqFeature, List[pysam.AlignedSegment]],
                             Dict[Tuple[SeqFeature.SeqFeature, SeqFeature.SeqFeature],
                                  List[pysam.AlignedSegment]]
                             ]:
    """
     * @description:
        Put reads (from uniq_bamis) to innerexon and junctions.
        1.  Reads are either defined to innerexon or junction
        2.  Reads map once, mapped sequence may be annotated to multipile exons
        3.  A read of junction can only annotated to pairs of reads:
            i.  in the same gene
            ii. have correct intron
     * @param {locats} generated by line_exons
     * @param {exon_loc} generated by line_exons
     * @param {anno_method} the method of annotation
     * @return {innerexon, junctions}
        innerexon: {exon: [valid reads, ]}
        junctions: {(left exon, right exon): [valid reads, ]}
    """
    def exonid_gene_split(anno_method):
        if anno_method == 'gencode':
            return ':'
        elif anno_method == 'refseq':
            return '-'

    split = exonid_gene_split(anno_method)
    innerexon: Dict[SeqFeature.SeqFeature, List[pysam.AlignedSegment]] = {}
    junctions: Dict[Tuple[SeqFeature.SeqFeature, SeqFeature.SeqFeature], List[pysam.AlignedSegment]] = {}
    failreads = set()
    for read_id, read in uniq_bamis.items():
        start = set(exon_iter(read.reference_start, locats, exon_loc))  # or position[0] - 1
        end = set(exon_iter(read.reference_end - 1, locats, exon_loc))  # or position[-1] - 1
        if not (start and end):
            failreads.add(read)
        l, m, r = start - end, start & end, end - start
        if 'N' in read.cigarstring:
            lend, rstart = read_junc(read)
            innergene: Dict[str, Tuple[SeqFeature.SeqFeature, SeqFeature.SeqFeature]] = {}
            for exon1 in l:
                if exon1.location.end.position < lend:
                    continue
                gene_name = exon1.id.split(split)[1]
                assert gene_name not in innergene
                innergene[gene_name] = exon1
            for exon2 in r:
                gene_name = exon2.id.split(split)[1]
                if gene_name not in innergene:
                    continue
                #assert not isinstance(innergene[gene_name], tuple)
                exon1: SeqFeature.SeqFeature = innergene[gene_name]
                if rstart < exon1.location.end.position:
                    continue
                innergene[gene_name] = exon1, exon2
            junction = {exons for exons in innergene.values() if isinstance(exons, tuple)}
            for exons in junction:
                junctions.setdefault(exons, []).append(read)
        else:
            for exon in m:
                innerexon.setdefault(exon, []).append(read)
    return innerexon, junctions, failreads


def analyze(fo=sys.stderr):
    with open() as gi:
        for seqfeature in GFF.parse(gi, base_dict={}):
            break
    locats, exon_loc = line_exons(extract_exons(seqfeature))

    exon_num = sum((len(exons) for exons in exon_loc.values()))
    print(f'total exon number: {exon_num}')

    uniq_bamis = load_bam_2_memory()
    innerexon, junctions, failreads = map_read_exon(
        uniq_bamis, locats, exon_loc, 'gencode')
    total_mapped_reads = len(uniq_bamis)  # 2818412
    print(f'Q2. '
          f'The number of reads that *can be aligned to* the human genome '
          f': {total_mapped_reads}')
    print('', total_mapped_reads, sep='\t', end='', file=fo)  # * report

    innerexon_reads = set()
    for reads in innerexon.values():
        innerexon_reads |= set(reads)
    print(f'Q4.a. '
          f'How many of the aligned reads *locate inside an exon* '
          f': {len(innerexon_reads)}')
    print('', len(innerexon_reads), sep='\t', end='', file=fo)  # * report

    junctions_reads = set()
    for reads in junctions.values():
        junctions_reads |= set(reads)
    print(f'Q4.b. '
          f'How many of the aligned reads *include exon junction* '
          f': {len(junctions_reads)}')
    print('', len(junctions_reads), sep='\t', end='', file=fo)  # * report


if __name__ == '__main__':
    report_file = 'Analyze/report.1.tsv'
    with open(report_file, 'w') as fo:
        analyze(fo)
