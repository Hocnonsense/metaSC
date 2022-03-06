# -*- coding: utf-8 -*-
"""
 * @Date: 2022-03-06 17:57:59
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-03-06 19:06:38
 * @FilePath: /metaSC/PyLib/test/PyLibTool/test_file_info.py
 * @Description:
    sth here
"""

from PyLib.PyLibTool.file_info import extract_doc, verbose_import, basicConfig


def test_extract_doc():
    doc_dict = extract_doc(__doc__)
    assert doc_dict["FilePath"] == "/metaSC/PyLib/test/PyLibTool/test_file_info.py"


def test_verbose_import():
    logger = verbose_import(__name__, __doc__)
    assert logger.name == __name__


def test_basicConfig():
    level = basicConfig()
    assert level == "INFO"
    level = basicConfig("DEBUG")
    assert level == basicConfig()
    assert level == "DEBUG"
