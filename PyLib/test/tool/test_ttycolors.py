# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-28 23:26:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 00:12:33
 * @FilePath: /metaSC/PyLib/test/tool/test_ttycolors.py
 * @Description:
"""
from PyLib.tool.ttycolors import color_text


def test_color_text():
    a = color_text("hellow color!", color="red", weight="bold")
    # assert a == "\x1b[0;31mhellow color!\x1b[0m"
    print(a)
