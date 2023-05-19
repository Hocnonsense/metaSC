# -*- coding: utf-8 -*-
"""
 * @Date: 2023-05-19 09:48:01
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-19 11:20:17
 * @FilePath: /metaSC/PyLib/test/biotool/test_download.py
 * @Description:
"""

import pytest
from PyLib.biotool import download
from PyLib.test import test_temp_path


@pytest.mark.slowdownload
def test_download_refseq():
    g = download.download_fna("GCF_000215995", test_temp_path)
    g1 = download.download_fna("GCF_000215995", test_temp_path)
    assert g.stat().st_mtime == g1.stat().st_mtime


@pytest.mark.slowdownload
def test_download_wgs():
    g = download.download_fna("BDQM00000000", test_temp_path)
    g1 = download.download_fna("BDQM00000000", test_temp_path)
    assert g.stat().st_mtime == g1.stat().st_mtime


@pytest.mark.slowdownload
def test_download_genebank():
    g = download.download_fna("CP002779", test_temp_path)
    g1 = download.download_fna("CP002779", test_temp_path)
    assert g.stat().st_mtime == g1.stat().st_mtime
