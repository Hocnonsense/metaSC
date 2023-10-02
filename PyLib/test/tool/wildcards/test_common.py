# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-30 21:39:00
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-30 21:44:46
 * @FilePath: /metaSC/PyLib/test/tool/wildcards/test_common.py
 * @Description:
"""

from PyLib.tool.wildcards import common as cmn


def test_common():
    frame = cmn.inspect.currentframe()
    assert cmn.stepout(0) == {}
    assert cmn.stepout(1)["frame"] == frame
    assert "frame" not in cmn.stepout(2)
