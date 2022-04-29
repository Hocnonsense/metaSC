# -*- coding: utf-8 -*-
"""
 * @Date: 2022-03-06 17:57:59
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 11:45:14
 * @FilePath: /metaSC/PyLib/test/PyLibTool/test_file_info.py
 * @Description:
    sth here
"""

from PyLib.PyLibTool.file_info import extract_doc, verbose_import, basicConfig, md5sum
from PyLib.test import test_file_path


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


def test_md5sum():
    file = test_file_path / "GCF_000215995.fna"
    assert (
        md5sum(
            text=file.open("rb").read(), file=file, c="a664eed1238f10bee8974c894f4b0191"
        )
        == "a664eed1238f10bee8974c894f4b0191"
    )
    assert not md5sum(bytes(file), c="a664eed1238f10bee8974c894f4b0191")
