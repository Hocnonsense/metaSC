# -*- coding: utf-8 -*-
"""
 * @Date: 2021-05-19 12:52:51
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-10 12:55:24
 * @FilePath: /metaSC/PyLib/reader/iters.py
 * @Description:
"""

from ast import literal_eval as eval
from typing import Iterable, List, Optional, Set, TextIO, Tuple, Union

from PyLib.biotool.kraken import kraken_level_filter
from PyLib.PyLibTool.file_info import verbose_import

verbose_import(__name__, __doc__)


def read_table(
    text: Union[List[str], str, TextIO],
    sep="\t",
    annot="#",
    title: Optional[List[str]] = None,
    openit=False,
):
    if openit:
        text = open(text)  # type: ignore
    for line in text:
        if line.startswith(annot):
            if title is not None:
                title.clear()
                title.extend(line[len(annot) :].rstrip().split(sep))
            continue
        values = line.strip().split(sep)
        if values:
            yield values
    if openit:
        text.close()  # type: ignore


def DASTool_summary_iter(text: TextIO) -> Iterable[List[str]]:
    """{'header': ['bin',
            'uniqueBacSCGs', 'redundantBacSCGs', 'uniqueArcSCGs', 'redundantArcSCGs',
            'bacRatio', 'arcRatio',
            'size', 'contigs', 'N50',
            'binScore', 'SCG_completeness', 'SCG_redundancy'
        ],
        'version': '2021-5-19 12:54:57'
    }"""
    header = eval(DASTool_summary_iter.__doc__)["header"]  # type: ignore
    for values in read_table(text):
        if values[0] == header[0]:
            assert all((header_i == value for header_i, value in zip(header, values)))
        if values:
            yield values


def DASTool_scaffolds2bin_iter(text: TextIO) -> Iterable[Tuple[str, Set[str]]]:
    """{
        'header': ['MAG_Id', 'Contigs']
    }"""
    MAG_Id = ""
    scaffolds: Set[str] = set()
    for scaffold, MAG in read_table(text):
        if MAG != MAG_Id:
            if MAG_Id:  # only work at very beginning
                yield MAG_Id, scaffolds
            MAG_Id = MAG
            scaffolds = set()
        scaffolds.add(scaffold)
    yield MAG_Id, scaffolds


def checkm_iter(text: TextIO) -> Iterable[Tuple[str, List[str]]]:
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
        if "\t" in line:
            break
        if line.startswith("--"):
            continue
        elif line.startswith("  Bin Id  "):
            header = line
            assert text.readline().startswith("--")
            line = text.readline()
            assert not line.startswith("  Bin Id  ")
        values = line.strip().split()
        yield values[0], [values[1], values[2][1:-1], values[3], *values[-3:]]
    else:
        return

    while line:
        if line.startswith("Bin Id"):
            header = line
            assert header.split()[12] in checkm_iter.__doc__  # type: ignore
            line = text.readline()
        values = line.strip().split()
        yield values[0], [values[1], values[2][1:-1], values[3], *values[-3:]]
        line = text.readline()


def gtdbtk_1_0_2_iter(text: TextIO):
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
    header = eval(gtdbtk_1_0_2_iter.__doc__)["in_header"]  # type: ignore
    for values in read_table(text):
        if values[0] == header[0]:
            assert all((header_i == value for header_i, value in zip(header, values)))
            continue
        yield values[0], [
            *values[1:4],
            values[7],
            values[12],
            values[15],
            values[17],
            values[18],
        ], [taxon.split("__")[1] for taxon in values[1].split(";")]


def gtdbtk_1_5_1_iter(text: TextIO):
    """{
        'in_header': [
            'user_genome',  # 0
            'classification', 'fastani_reference', 'fastani_reference_radius',  # 3
            'fastani_taxonomy', 'fastani_ani', 'fastani_af',  # 6
            'closest_placement_reference',  # 7
            'closest_placement_radius', 'closest_placement_taxonomy', 'closest_placement_ani', 'closest_placement_af',  # 11
            'pplacer_taxonomy', 'classification_method', 'note', 'other_related_references(genome_id,species_name,radius,ANI,AF)',  # 15
            'msa_percent', 'translation_table', 'red_value', 'warnings'  # 19
        ],
        'out_header': (
            ['user_genome'],
            [
                'classification', 'fastani_reference', 'fastani_reference_radius',
                'closest_placement_reference', 'classification_method', 'msa_percent',
                'red_value', 'warnings'
            ],
            ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        )
    }"""
    header = eval(gtdbtk_1_5_1_iter.__doc__)["in_header"]  # type: ignore
    for values in read_table(text):
        if values[0] == header[0]:
            assert all((header_i == value for header_i, value in zip(header, values)))
            continue
        yield values[0], [
            *values[1:4],
            values[7],
            values[13],
            values[16],
            values[18],
            values[19],
        ], [taxon.split("__")[1] for taxon in values[1].split(";")]


gtdbtk_iter = gtdbtk_1_5_1_iter


def featureCounts_iter(text: TextIO, ititle: Optional[List] = None):
    """ read <in_file> of featureCounts by:
            featureCounts \\
                -a ${gff} \\
                -o <in_file> \\
                -t CDS -g ID \\
                -p ${bam}
       @return:
           |   0   | 1  |  2   | 3  |   4   |   5   |  6 | ...
            Geneid, Chr, Start, End, Strand, Length, bam1  ...
    """
    for values in read_table(text):
        # drop title line
        if values[0].startswith("Geneid") and isinstance(ititle, list):
            ititle.clear()
            ititle.extend((*values[:6], values[6:]))
        Geneid, Chr, Start, End, Strand, Length = values[:6]
        bam = values[6:]
        yield (
            Geneid,
            Chr,
            int(Start),
            int(End),
            Strand,
            int(Length),
            [int(bami) for bami in bam],
        )


def emapper_g5f52332_iter(text: TextIO):
    """{
        'header': [
            'query_name', 'seed_eggNOG_ortholog', 'seed_ortholog_value', 'seed_ortholog_score',
            'Predicted_taxonomic_group', 'Predicted_protein_name', 'Gene_Ontology_terms',
            'EC_number', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC',
            'CAZy', 'BiGG_Reaction', 'tax_scope',
            'eggNOG_OGs', 'bestOG', 'COG_Functional_Category', 'eggNOG_free_text'
        ],
        'version': 'emapper-2.1.3-ee27b8e'
    }"""
    for values in read_table(text):
        yield values


def emapper_ee27b8e_iter(text: TextIO):
    """{
        'header': [
            'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
            'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
            'GOs', 'EC',
            'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC',
            'CAZy', 'BiGG_Reaction', 'PFAMs'
        ],
        'version': 'emapper-2.1.3-ee27b8e'
    }"""
    for values in read_table(text):
        yield values


emapper_iter = emapper_ee27b8e_iter


kraken_iter = kraken_level_filter
