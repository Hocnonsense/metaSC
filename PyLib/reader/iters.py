# -*- coding: utf-8 -*-
"""
 * @Date: 2021-05-19 12:52:51
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-06-07 17:40:17
 * @FilePath: /2021_05-MT10kSW/Scripts/PyLib/reader/iters.py
 * @Description:
"""

from ast import literal_eval as eval
from io import FileIO
from typing import Iterable, List, Tuple

from PyLib.PyLibTool.file_info import verbose_import


verbose_import(__name__, __doc__)


def DASTool_summary_iter(text: FileIO):
    """{'header': ['bin',
            'uniqueBacSCGs', 'redundantBacSCGs', 'uniqueArcSCGs', 'redundantArcSCGs',
            'bacRatio', 'arcRatio',
            'size', 'contigs', 'N50',
            'binScore', 'SCG_completeness', 'SCG_redundancy'
        ],
        'version': '2021-5-19 12:54:57'
    }"""
    header = text.readline()
    assert header == '\t'.join(eval(DASTool_summary_iter.__doc__)['header'])
    for line in text:
        values = line.strip().split('\t')
        if values:
            yield values


def checkm_iter(text: FileIO) -> Iterable[Tuple[str, List[str]]]:
    """{
        'header': (
            'Bin Id',
            [
                'Marker lineage', '# genomes', '# markers', '# marker sets',
                '0', '1', '2', '3', '4', '5+',
                'Completeness', 'Contamination', 'Strain heterogeneity'
            ]
        )
    }"""
    for line in text:
        if '\t' in line:
            break
        if line.startswith('--'):
            continue
        elif line.startswith('  Bin Id  '):
            header = line
            assert text.readline().startswith('--')
            line = text.readline()
            assert not line.startswith('  Bin Id  ')
        values = line.strip().split()
        yield values[0], [
            values[1], values[2][1:-1], values[3],
            *values[-3:]
        ]
    else:
        return

    while line:
        if line.startswith('Bin Id'):
            header = line
            assert header.split()[12] in checkm_iter.__doc__
            line = text.readline()
        values = line.strip().split()
        yield values[0], [
            values[1], values[2][1:-1], values[3],
            *values[-3:]
        ]
        line = text.readline()


def gtdbtk_iter(text: FileIO):
    """{
        'in_header': [
            'user_genome',  # 0
            'classification', 'fastani_reference', 'fastani_reference_radius',  # 3
            'fastani_taxonomy', 'fastani_ani', 'fastani_af',  # 6
            'closest_placement_reference',  #7
            'closest_placement_taxonomy', 'closest_placement_ani', 'closest_placement_af',  #10
            'pplacer_taxonomy', 'classification_method', 'note', 'other_related_references(genome_id,species_name,radius,ANI,AF)',
            'aa_percent', 'translation_table', 'red_value', 'warnings'
        ],
        'out_header': (
            ['user_genome'],
            [
                'classification', 'fastani_reference', 'fastani_reference_radius',
                'closest_placement_reference', 'classification_method', 'aa_percent',
                'red_value', 'warnings'
            ],
            ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        )
    }"""
    header = text.readline()
    assert header == '\t'.join(eval(gtdbtk_iter.__doc__)['in_header']) + '\n'  # 12
    for line in text:
        if not line.strip():
            continue
        values = line.strip().split('\t')
        yield values[0], [
            *values[1:4],
            values[7], values[12], values[15],
            values[17], values[18]
        ], [
            taxon.split('__')[1]
            for taxon in values[1].split(';')
        ]
