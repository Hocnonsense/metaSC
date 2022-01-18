# -*- coding: utf-8 -*-
"""
 * @Date: 2021-07-31 15:28:13
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-18 11:12:42
 * @FilePath: /metaSC/PyLib/biotool/blast.py
 * @Description:
"""

import xml.dom.minidom
from io import FileIO
from typing import List, Tuple

from Bio import Seq

BLAST_XML_Hit_TITLE = [
    'Hit_def', 'Hit_accession', 'Hit_len'
]

BLAST_XML_Hsp_TITLE = [
    'Hsp_num', 'Hsp_bit-score', 'Hsp_score', 'Hsp_evalue',
    'Hsp_query-from', 'Hsp_query-to', 'Hsp_hit-from', 'Hsp_hit-to',
    'Hsp_query-frame', 'Hsp_hit-frame',
    'Hsp_identity', 'Hsp_positive', 'Hsp_gaps', 'Hsp_align-len',
    'Hsp_qseq', 'Hsp_hseq', 'Hsp_midline'
]


def getEle(xmlEle: xml.dom.minidom.Element, TagName):
    Element: xml.dom.minidom.Element = xmlEle.getElementsByTagName(TagName)[0]
    firstChild: xml.dom.minidom.Text = Element.firstChild
    return firstChild.data


def blast_xml_iter(blast_file: FileIO):
    '''{
        title: ['query_def',  # 0
                # BLAST_XML_Hit_TITLE
                'Hit_def', 'Hit_accession', 'Hit_len',  # 3
                # BLAST_XML_Hsp_TITLE
                'Hsp_num', 'Hsp_bit-score', 'Hsp_score', 'Hsp_evalue',  # 7
                'Hsp_query-from', 'Hsp_query-to', 'Hsp_hit-from', 'Hsp_hit-to',  # 11
                'Hsp_query-frame', 'Hsp_hit-frame',  # 13
                'Hsp_identity', 'Hsp_positive', 'Hsp_gaps', 'Hsp_align-len',  # 17
                'Hsp_qseq', 'Hsp_hseq', 'Hsp_midline'
        ]
    }'''
    DOMTree: xml.dom.minidom.Document = xml.dom.minidom.parse(blast_file)

    blast_results: xml.dom.minidom.Element = DOMTree.documentElement
    BlastOutput_iterations: xml.dom.minidom.Element = blast_results.getElementsByTagName(
        'BlastOutput_iterations')[0]

    for Iteration in BlastOutput_iterations.getElementsByTagName('Iteration'):
        query_def = getEle(Iteration, 'Iteration_query-def')
        for Hit in Iteration.getElementsByTagName('Iteration_hits')[0].getElementsByTagName('Hit'):
            Hit_values = [getEle(Hit, key) for key in BLAST_XML_Hit_TITLE]
            for Hsp in Hit.getElementsByTagName('Hit_hsps')[0].getElementsByTagName('Hsp'):
                Hsp_values = [getEle(Hsp, key) for key in BLAST_XML_Hsp_TITLE]
                yield [query_def] + Hit_values + Hsp_values


def blast_unmatch_iter(values: List[str], wind=10) -> List[Tuple[List[str], List[int], List[int]]]:
    """ read one of blast alignment and output each mismatch
    * @param {List[str]} values: values itered from blast_iter
    * @param {int} wind: length of wind sequence
    * @return List[str]: alignment of hit seq
        [
            mismatch, base_from, base_to:
                mismatch: enum(insert, delete, mismatch) compared to hit:
                    insert means query sequence have one more base at the unmatched position
            hit_from, hit_at, hit_to:
                absolute 0-start index of hit sequence;
            hsp_qseq, hsp_midl, hsp_hseq:
                alignment of unmatched segment;
        ]
    """
    wind, half_wind = int(wind), int(wind) // 2

    def lwind(index, hseqstart):
        return (index - wind + half_wind + hseqstart) // wind * wind - hseqstart + 1

    def rwind(index, hseqstart):
        return (index + wind + half_wind + hseqstart) // wind * wind - hseqstart + 1

    def extend_unmatch(index, newindex=0):
        index, newindex = min(index, newindex), max(index, newindex)
        newunmatch = Hsp_hseq.count('-', index, newindex)
        while newunmatch:
            index, newindex = newindex, newindex + newunmatch
            newunmatch = Hsp_hseq.count('-', index, newindex)
        return newindex

    Hsp_hit_from: int = int(values[10])
    Hsp_hit_to: int = int(values[11])

    Hsp_qseq, Hsp_hseq, Hsp_midline = values[-3:]
    if Hsp_hit_to < Hsp_hit_from:
        Hsp_hit_from, Hsp_hit_to = Hsp_hit_to, Hsp_hit_from
        Hsp_qseq = Seq.Seq(Hsp_qseq).reverse_complement()
        Hsp_hseq = Seq.Seq(Hsp_hseq).reverse_complement()
        Hsp_midline = ''.join(reversed(Hsp_midline))
    Hsp_hseq_len = len(Hsp_hseq.replace('-', ''))

    unmatch = Hsp_midline.find(' ', 0)
    while unmatch != -1:
        if Hsp_qseq[unmatch] == '-':
            mismatch_type = 'delete'
        elif Hsp_hseq[unmatch] == '-':
            mismatch_type = 'insert'
        else:
            mismatch_type = 'mismatch'

        ref_unmatch = unmatch - Hsp_hseq.count('-', 0, unmatch)
        ref_lextend = max(lwind(ref_unmatch, Hsp_hit_from), 0)
        ref_rextend = min(rwind(ref_unmatch, Hsp_hit_from), Hsp_hseq_len)

        lextend = extend_unmatch(ref_lextend)
        rextend = extend_unmatch(ref_rextend)

        assert len(Hsp_hseq[lextend: rextend].replace('-', '')) == ref_rextend - ref_lextend

        yield ([mismatch_type, Hsp_hseq[unmatch], Hsp_qseq[unmatch]],
               [Hsp_hit_from + index for index in
                (ref_lextend, ref_unmatch, ref_rextend - 1)],
               [seq[lextend: rextend] for seq in
                (Hsp_qseq, Hsp_midline, Hsp_hseq)])

        unmatch = Hsp_midline.find(' ', unmatch + 1)
