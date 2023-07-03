# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-24 10:50:36
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-07-03 14:01:20
 * @FilePath: /metaSC/PyLib/setup.py
 * @Description:
    setup
"""

import os
from pathlib import Path

os.chdir(Path(__file__).parent.parent)


from setuptools import setup, find_packages

setup(
    name="metaSC",
    author="Hwrn",
    author_email="Hwrn.aou@sjtu.edu.cn",
    packages=find_packages(),
    version="0.0.3",
)
