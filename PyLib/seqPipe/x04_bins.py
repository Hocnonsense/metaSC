# -*- coding: utf-8 -*-
"""
 * @Date: 2021-08-14 14:35:18
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-08-14 17:22:47
 * @FilePath: /metaSC/PyLib/seqPipe/x04_bins.py
 * @Description:
"""
import os
from typing import Dict, List, Set, TextIO, Tuple

from PyLib.reader.read_outputs import fasta, jgi_depths
from PyLib.reader.iters import checkm_iter, gtdbtk_iter
from PyLib.biotool.fna_msg import statistic_fna, seq_total_depth
from PyLib.PyLibTool.file_info import verbose_import

logger = verbose_import(__name__, __doc__)


def MAG_seq_features(bin_filename: str,
                     ctg_depth: Dict[str, Tuple[Tuple[int, float],
                                                List[float], List[float]]] = []):
    """{
        'values': [
            'SeqNumbers', 'MaxLength', 'GenomeSize',
            'GC', 'N50', 'L50',
            'total_depth', '{depths}'
        ],
        'keys': {
            'ctg_depth': 'generated by jgi_depths from PyLib.reader.read_outputs'
        }
    }"""
    fna_msg = statistic_fna(fasta(bin_filename))
    totalAvgDepth, depths = seq_total_depth(ctg_depth, fasta(bin_filename))
    return list(fna_msg) + [totalAvgDepth] + depths


def wrap_depth(fi: TextIO):
    title = next(fi).split()
    yield '\t'.join(title + [title[-1] + '-var'])
    for line in fi:
        yield line.replace('\n', '\t0\n')


def merge_checkm_bin(basedir: str, metawrap: bool = False,
                     summary_dict: Dict[str, List[str]] = None):
    """{
        'key': 'Bin_Id',
        'values': [
            'Marker_lineage', 'lineage_UID', 'genomes',  # checkm
            'Completeness', 'Contamination', 'heterogeneity',
            'SeqNumbers', 'MaxLength', 'GenomeSize', 'GC', 'N50', 'L50',
            'totalAvgDepth', '{depths}',
        ],
        'values+': [
            'SeqNumbers', 'MaxLength', 'GenomeSize', 'GC', 'N50', 'L50',
            'totalAvgDepth', '{depths}',
        ]
    }"""
    DEPTH_FILE_PATH = f'{basedir}/work_files/mb2_master_depth.txt'
    if metawrap:
        FORMAT_BIN_FILE_PATH = '_bins'
        with open(DEPTH_FILE_PATH) as fi:
            sample_list, ctg_depth = jgi_depths(wrap_depth(fi))
    else:
        FORMAT_BIN_FILE_PATH = '_DASTool_bins'
        with open(DEPTH_FILE_PATH) as fi:
            sample_list, ctg_depth = jgi_depths(fi)

    # read checkm
    if summary_dict is None:
        file_checkm = f'{basedir}/checkms/report.tsv'
        with open(file_checkm) as fi:
            for MAG, values in checkm_iter(fi):
                summary_dict[MAG] = values

    for binId, values in summary_dict.items():
        genome_features = MAG_seq_features(f'{basedir}/{FORMAT_BIN_FILE_PATH}/{binId}.fa', ctg_depth)
        summary_dict[binId] = values + genome_features
    return summary_dict


def merge_checkm_taxon(basedir: str, file_taxon: str,
                       summary_dict: Dict[str, List[str]] = None):
    """{
        'key': 'Bin_Id',
        'values': [
            'Marker_lineage', 'lineage_UID', 'genomes',  # checkm
            'Completeness', 'Contamination', 'heterogeneity',
            'classification', 'fastani_reference', 'fastani_reference_radius',
            'closest_placement_reference', 'classification_method', 'aa_percent',
            'red_value', 'warnings',
            'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'
        ],
        'values+': [
            'classification', 'fastani_reference', 'fastani_reference_radius',
            'closest_placement_reference', 'classification_method', 'aa_percent',
            'red_value', 'warnings',
            'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'
        ]
    }"""
    if summary_dict is None:
        file_checkm = f'{basedir}/checkms/report.tsv'
        with open(file_checkm) as fi:
            for MAG, values in checkm_iter(fi):
                summary_dict[MAG] = values

    for kindom in ['ar122', 'bac120']:
        filename = f'{file_taxon}/gtdbtk.{kindom}.summary.tsv'
        if not os.path.exists(filename):
            logger.warning(f'{filename} not exists')
        with open(filename) as fi:
            for MAG, values, taxon in gtdbtk_iter(fi):
                summary_dict[MAG] += values + taxon
    return summary_dict
