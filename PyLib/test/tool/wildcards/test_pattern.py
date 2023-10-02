# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-26 11:55:48
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 16:15:15
 * @FilePath: /metaSC/PyLib/test/tool/wildcards/test_pattern.py
 * @Description:
"""

import re
from pathlib import Path

import pytest

from PyLib.tool.wildcards import pattern as ptn

test_patterns = [
    ptn.PathPattern(i)  # type: ignore
    for i in (
        "README.md",
        "pipe/{site}/02_assem_{site}.fa",
        Path("pipe/{site}/02_mapping_{site}..{layer}.fa"),
        "pipe/{site}/03_annot_{site}-{annot,kofam|eggnog|mantis}.faa",
    )
]


@pytest.mark.xfail
def test_create_pattern_blank():
    ptn.PathPattern("{sth} like {}")


@pytest.mark.xfail
def test_create_pattern_number():
    ptn.PathPattern("{sth} like {1}")


def test_get_wildcard_names():
    assert test_patterns[0].names == []
    assert test_patterns[1].names == ["site"]
    assert test_patterns[2].names == ["site", "layer"]
    assert test_patterns[3].names == ["site", "annot"]


def test_get_wildcards_complex():
    p4 = ptn.PathPattern("a: {a[0]}| b: {b}")
    assert p4.names == ["b"]


def test_get_wildcard_constraints():
    assert not test_patterns[0].constraints
    assert not test_patterns[1].constraints
    assert not test_patterns[2].constraints
    assert test_patterns[3].constraints == {"annot": re.compile("kofam|eggnog|mantis")}


def test_strip_wildcard_constraints():
    assert test_patterns[0].strip_constraints == test_patterns[0].pattern
    assert test_patterns[1].strip_constraints == test_patterns[1].pattern
    assert test_patterns[2].strip_constraints == Path(test_patterns[2].pattern)
    assert (
        test_patterns[3].strip_constraints == "pipe/{site}/03_annot_{site}-{annot}.faa"
    )


def test_constraint_replace():
    assert test_patterns[1].replace_constraint({"site": "M*_2"}).constraints == {
        "site": re.compile("M*_2")
    }
    assert test_patterns[1].replace_constraint({"layer": "M*_2"}).constraints == {}
    assert test_patterns[2].replace_constraint({"layer": "M*_2"}).constraints == {
        "layer": re.compile("M*_2")
    }


def test_regex():
    assert test_patterns[0].regex == re.compile("README\\.md$")
    assert test_patterns[1].regex == re.compile(
        "pipe/(?P<site>.+)/02_assem_(?P=site)\\.fa$"
    )
    assert test_patterns[2].regex == re.compile(
        "pipe/(?P<site>.+)/02_mapping_(?P=site)\\.\\.(?P<layer>.+)\\.fa$"
    )
    assert test_patterns[3].regex == re.compile(
        "pipe/(?P<site>.+)/03_annot_(?P=site)\\-(?P<annot>kofam|eggnog|mantis)\\.faa$"
    )


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


def test_glob_wildcards():
    assert not list(test_patterns[0].glob(files))
    wc = list(test_patterns[1].glob(files))
    assert [i[1].site for i in wc] == ["M1", "M2"]
    wc = list(test_patterns[2].glob(files))
    assert wc[0][1].site == "M1" and wc[0][1].layer == "0-2" and len(wc) == 1
    wc = list(test_patterns[3].glob(files))
    assert wc[0][1].site == "M1" and wc[0][1].annot == "mantis" and len(wc) == 1
