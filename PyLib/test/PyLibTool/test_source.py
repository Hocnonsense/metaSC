# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-04 20:02:28
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-04 20:46:41
 * @FilePath: /metaSC/PyLib/test/PyLibTool/test_source.py
 * @Description:
"""

import os
from PyLib.PyLibTool.source import source, source_env

a = 1


def test_source():
    assert source(__file__).a == 1


def test_source_env():
    assert os.system(f"test={__file__} python {__file__}") == 0
    assert os.system(f"test='11' python {__file__}") != 0


if __name__ == "__main__":
    print(source_env("test").a)
    a = 2
