# -*- coding: utf-8 -*-
"""
 * @Date: 2022-03-06 16:51:44
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-03-06 19:11:24
 * @FilePath: /metaSC/PyLib/test/PyLibTool/test_tmpPkl.py
 * @Description:
"""

import os
from PyLib.PyLibTool.tmpPkl import TmpPkl


file_path = os.path.join(os.path.dirname(__file__), "test.pkl")
if os.path.exists(file_path):
    os.remove(file_path)


@TmpPkl(file_path, force_rewrite=False)
def atest(i, j):
    """Test i+j."""
    return i + j


def test_funct():
    a = atest(1, 2)
    assert a == 3
    a = atest(2, 2)
    assert a == 3


def test_with():
    with TmpPkl(
        "test.pkl", desc="with i + j", force_rewrite=True, situ=__file__
    ) as tmp1:
        tmp1.desc = "with yybsy"
        if tmp1.force_rewrite:
            tmp1.last_results = "ii" + "Vi"
            tmp1.desc = {"i": "ii", "j": "vi"}
            word = tmp1.last_results
        assert tmp1.meta["desc"] == {"i": "ii", "j": "vi"}
        assert tmp1.cache[atest] == True
    assert word == "iiVi"


def test_delete():
    if os.path.exists(file_path):
        os.remove(file_path)