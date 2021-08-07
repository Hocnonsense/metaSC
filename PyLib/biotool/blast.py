# -*- coding: utf-8 -*-
"""
 * @Date: 2021-07-31 15:28:13
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-08-07 22:27:34
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


def blast_unmatch_iter(values: List[str], wind=10, foreach=False) -> List[Tuple[str, List[int], List[int]]]:
    """ read one of blast alignment and output each mismatch
    * @param {List[str]} values: values itered from blast_iter
    * @param {int} wind: length of wind sequence
    * @param {bool} foreach:
        ? WRANING: not implemented yet.
        if True: return each alignment;
        else: if mismatch is closed to another, they will be itered togather
    * @return List[str]: alignment of hit seq
        [
            mismatch: enum(insert, delete, mismatch) compared to hit:
                insert means query sequence have one more base at the unmatched position
            hit_from, hit_at, hit_to:
                absolute 0-start index of hit sequence;
            hsp_qseq, hsp_midl, hsp_hseq:
                alignment of unmatched segment;
        ]
    """
    Hsp_hit_from: int = int(values[10])
    Hsp_hit_to: int = int(values[11])
    Hsp_align_len = int(values[17])
    Hsp_qseq, Hsp_hseq, Hsp_midline = values[-3:]
    if Hsp_hit_to < Hsp_hit_from:
        Hsp_hit_from, Hsp_hit_to = Hsp_hit_to, Hsp_hit_from
        Hsp_qseq = Seq.Seq(Hsp_qseq).reverse_complement()
        Hsp_hseq = Seq.Seq(Hsp_hseq).reverse_complement()
        Hsp_midline = ''.join(reversed(Hsp_midline))
    length = Hsp_midline.find(' ', 0)
    wind, half_wind = int(wind), wind // 2 + 1
    while length != -1:
        extendleft = max((length - wind + half_wind) // wind * wind + 1,
                         0)
        extendright = min((length + wind + half_wind) // wind * wind + 1,
                          Hsp_align_len)
        if Hsp_qseq[length] == '-':
            mismatch_type = 'delete'
        elif Hsp_hseq[length] == '-':
            mismatch_type = 'insert'
        else:
            mismatch_type = 'mismatch'
        yield (
            mismatch_type,
            [Hsp_hit_from + index - 1 - Hsp_hseq.count('-', 0, index)
             for index in (extendleft, length + 1, extendright)],
            [seq[extendleft: extendright] for seq in
             (Hsp_qseq, Hsp_midline, Hsp_hseq)]
        )
        length = Hsp_midline.find(' ', length + 1)
