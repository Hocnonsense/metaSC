# -*- coding: utf-8 -*-
"""
 * @Date: 2021-05-03 20:18:57
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-05-06 20:11:55
 * @FilePath: /metaSC/PyLib/reader/DASTool_summary.py
 * @Description:
"""

from ast import literal_eval as eval
from io import FileIO


def summaryiter(text: FileIO):
    """['bin',
        'uniqueBacSCGs', 'redundantBacSCGs', 'uniqueArcSCGs', 'redundantArcSCGs',
        'bacRatio', 'arcRatio',
        'size', 'contigs', 'N50',
        'binScore', 'SCG_completeness', 'SCG_redundancy'
    ]"""
    header = text.readline()
    assert header == '\t'.join(eval(summaryiter.__doc__))
    for line in text:
        values = line.strip().split('\t')
        if values:
            yield values
