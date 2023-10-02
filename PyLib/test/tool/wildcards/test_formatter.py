# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-25 21:35:40
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 18:51:33
 * @FilePath: /metaSC/PyLib/test/tool/wildcards/test_formatter.py
 * @Description:
"""

from pathlib import Path
from PyLib.tool.wildcards import formatter as fmt
from PyLib.tool.wildcards.common import Pattern


a = "a1"
b = "b1"


def test_formatter_args():
    assert fmt.StrPattern("a: {}| b: {b}").xformat()("a 3") == "a: a 3| b: b1"


def test_formatter_kwargs():
    assert fmt.StrPattern("a: {a}| b: {b}").xformat()(a="a 3") == "a: a 3| b: b1"


def test_formatter_stepout():
    a = "a2"
    assert fmt.StrPattern("a: {a}| b: {b}").xformat()() == "a: a2| b: b1"

    def step1(stepout: int):
        a = "a4"

        def step2():
            a = "a5"
            return fmt.StrPattern("a: {a}| b: {b}").xformat(nstep=stepout)

        return step2()

    assert step1(0)(a="a0", b="b0") == "a: a0| b: b0"
    assert step1(1)() == "a: a5| b: b1"
    assert step1(2)() == "a: a4| b: b1"
    assert step1(3)() == "a: a2| b: b1"


def test_formatter_missing():
    assert (
        fmt.StrPattern("a: {a}| c: {c}").xformat(keep_missing=True)() == "a: a1| c: {c}"
    )


def test_formatter_quoter():
    assert fmt.StrPattern("a: {a}| b: {b}").xformat()(a=[1, 2, 3]) == "a: 1 2 3| b: b1"
    assert (
        fmt.StrPattern("a: {a}| b: {b}").xformat()(a={1: 2, 3: 4})
        == "a: 1=2,3=4| b: b1"
    )
    assert (
        fmt.StrPattern("a: {a}| b: {b}").xformat()(a={1: 2, 3: 4}.keys())
        == "a: dict_keys([1, 3])| b: b1"
    )
    assert (
        fmt.StrPattern("a: {a}| b: {b}").xformat(quote_all=True)(a="1 2 3")
        == "a: '1 2 3'| b: b1"
    )
    assert (
        fmt.StrPattern("a: {a}| b: {b}").xformat(quote_all=True)(a=[1, 2, 3])
        == "a: 1 2 3| b: b1"
    )


def test_expand():
    assert fmt.StrPattern("a: {a}| b: {b}").xexpand()(a="a1", b=b) == ["a: a1| b: b1"]
    assert fmt.StrPattern("a: {a}| b: {b}").xexpand()(a=["a1", "a2"], b=b) == [
        "a: a1| b: b1",
        "a: a2| b: b1",
    ]
    assert fmt.StrPattern("a: {a}| b: {b}").xexpand(keep_missing=True)(
        a=["a1", "a2"]
    ) == [
        "a: a1| b: {b}",
        "a: a2| b: {b}",
    ]
    assert fmt.StrPattern("a: {a}| b: {b}").xexpand(keep_missing=True)(
        a=["a1", ["a21", "a22"]]
    ) == [
        "a: a1| b: {b}",
        "a: ['a21', 'a22']| b: {b}",
    ]


def test_expand_args():
    assert fmt.StrPattern("a: {}| b: {b}").xexpand()(["a1", "a2"], b=b) == [
        "a: a1| b: b1",
        "a: a2| b: b1",
    ]
    assert fmt.StrPattern("a: {1}| b: {b}").xexpand()(["a0"], ["a1", "a2"], b=b) == [
        "a: a1| b: b1",
        "a: a2| b: b1",
    ]


def test_formatter_kwargs_complex():
    assert fmt.StrPattern("a: {a[0]}| b: {b}").xformat()(a=["a1"]) == "a: a1| b: b1"
    assert fmt.StrPattern("a: {a}| b: {b}").xexpand()(a=["a1"], b=b) == ["a: a1| b: b1"]
    assert fmt.StrPattern("a: {a[0]}| b: {b}").xexpand()(a=[["a1"], "a2"], b=b) == [
        "a: a1| b: b1",
        "a: a| b: b1",
    ]


def test_apply_wildcards():
    wcds = {"site": "M1", "layer": "0-2"}
    from PyLib.test.tool.wildcards.test_pattern import test_patterns

    assert fmt.StrPattern(test_patterns[0]).xformat(nstep=0)(**wcds) == "README.md"
    assert (
        fmt.StrPattern(test_patterns[1]).xformat(nstep=0)(**wcds)
        == "pipe/M1/02_assem_M1.fa"
    )
    assert isinstance(test_patterns[2], Pattern)
    assert isinstance(test_patterns[2]._type(), Path)
    assert isinstance(fmt.StrPattern(test_patterns[2])._type(), Path)
    assert fmt.StrPattern(test_patterns[2]).xformat(nstep=0)(**wcds) == test_patterns[
        2
    ]._type("pipe/M1/02_mapping_M1..0-2.fa")
    assert (
        fmt.StrPattern(test_patterns[3].strip_constraints).xformat()(
            **wcds, annot="mantis"
        )
        == "pipe/M1/03_annot_M1-mantis.faa"
    )
    assert (
        fmt.StrPattern(test_patterns[3].strip_constraints).xformat()(
            **wcds, annot="prodigal"
        )
        == "pipe/M1/03_annot_M1-prodigal.faa"
    )
    assert fmt.StrPattern(test_patterns[3].strip_constraints).xformat(
        keep_missing=True
    )(**wcds) == test_patterns[3]._type("pipe/M1/03_annot_M1-{annot}.faa")
    assert (
        fmt.StrPattern(test_patterns[3].strip_constraints).xformat(fill_missing=True)(
            **wcds
        )
        == "pipe/M1/03_annot_M1-.faa"
    )
    assert (
        fmt.StrPattern(test_patterns[3].strip_constraints).xformat(
            fill_missing=True,
            missing_value="{any}",
        )(**wcds)
        == "pipe/M1/03_annot_M1-{any}.faa"
    )
