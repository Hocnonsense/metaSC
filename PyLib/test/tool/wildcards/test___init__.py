# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-26 11:55:48
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 16:15:44
 * @FilePath: /metaSC/PyLib/test/tool/wildcards/test___init__.py
 * @Description:
"""

from pathlib import Path

import pytest

from PyLib.tool.wildcards import fmter, epder, glob_path, strip_constraints

test_patterns = [
    "README.md",
    "pipe/{site}/02_assem_{site}.fa",
    Path("pipe/{site}/02_mapping_{site}..{layer}.fa"),
    "pipe/{site}/03_annot_{site}-{annot,kofam|eggnog|mantis}.faa",
]

files = [
    "README.md",
    "pipe/M1/02_assem_M1.fa",
    "pipe/M1/02_assem_M2.fa",
    "pipe/M2/02_assem_M2.fa",
    "pipe/M1/02_mapping_M1..0-2.fa",
    "pipe/M1/03_annot_M1-mantis.faa",
    "pipe/M1/03_annot_M1-prodigal.faa",
    "pipe/M1/03_annot_M1-{annot}.faa",
]


def test_glob_path():
    assert not list(glob_path(test_patterns[0], files))
    wc = list(glob_path(test_patterns[1], files))
    assert [i[1].site for i in wc] == ["M1", "M2"]
    wc = list(glob_path(test_patterns[2], files))
    assert wc[0][1].site == "M1" and wc[0][1].layer == "0-2" and len(wc) == 1
    wc = list(glob_path(test_patterns[3], files))
    assert wc[0][1].site == "M1" and wc[0][1].annot == "mantis" and len(wc) == 1


def test_fmter():
    assert fmter(test_patterns[0])() == "README.md"
    assert (
        fmter(test_patterns[1], keep_missing=True)() == "pipe/{site}/02_assem_{site}.fa"
    )
    site = "TY"
    assert fmter(test_patterns[2], nstep=0)(site="FDZ", layer="0-2") == Path(
        "pipe/FDZ/02_mapping_FDZ..0-2.fa"
    )
    assert (
        fmter(strip_constraints(test_patterns[3]), nstep=1)(annot="ncbi")
        == "pipe/TY/03_annot_TY-ncbi.faa"
    )


@pytest.mark.xfail
def test_fmter_fail_on_pattern():
    fmter(test_patterns[3], nstep=1)(annot="ncbi")
