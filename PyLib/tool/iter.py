# -*- coding: utf-8 -*-
"""
 * @Date: 2022-02-08 22:20:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-08 22:20:36
 * @FilePath: /metaSC/PyLib/tool/iter.py
 * @Description:
"""


def better_grouper(inputs, n):
    iters = [iter(inputs)] * n
    return zip(*iters)
