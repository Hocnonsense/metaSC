# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-28 21:31:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-28 21:59:00
 * @FilePath: /metaSC/PyLib/test/tool/test_logger.py
 * @Description:
"""

import os
import sys

import PyLib.tool.logger as mylogger
from PyLib.test import test_file_path


def test_logger():
    with mylogger.Redirector():
        print(1)

    with mylogger.Redirector(
        stderr=open(test_file_path / "logger.err", "w"),
        stdout=open(test_file_path / "logger.out", "w"),
    ):
        print(2)
        os.system("echo 11 >&2")


def test_tee():
    sys.stdout = mylogger.Tee("w", filename=test_file_path / "tee.log")
    print(3)
