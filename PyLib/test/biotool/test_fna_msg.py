# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-29 11:58:01
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 12:03:17
 * @FilePath: /metaSC/PyLib/test/biotool/test_fna_msg.py
 * @Description:
"""

from PyLib.biotool import fna_msg
from PyLib.reader.read_outputs import fasta
from PyLib.test import test_file_path


def test_statistic_fna():
    genome = test_file_path / "GCF_000215995.fna"
    assert fna_msg.statistic_fna(fasta(genome)) == (
        1,
        1716818,
        1716818,
        0.5163750613052752,
        1716818,
        1,
    )
