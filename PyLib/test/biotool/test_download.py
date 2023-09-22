# -*- coding: utf-8 -*-
"""
 * @Date: 2023-05-19 09:48:01
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-01 12:31:46
 * @FilePath: /metaSC/PyLib/test/biotool/test_download.py
 * @Description:
"""

from pathlib import Path
import pytest
from PyLib.biotool import download
from PyLib.test import test_temp_path


@pytest.mark.slowdownload
def test_download_refseq():
    g = download.download_genome("GCF_000215995", test_temp_path)
    g1 = download.download_genome("GCF_000215995", test_temp_path)
    assert g.stat().st_mtime == g1.stat().st_mtime


@pytest.mark.slowdownload
def test_download_wgs():
    g = download.download_genome("BDQM00000000", test_temp_path)
    g1 = download.download_genome("BDQM00000000", test_temp_path)
    assert g.stat().st_mtime == g1.stat().st_mtime


@pytest.mark.slowdownload
def test_download_genebank():
    g = download.download_genome("CP002779", test_temp_path)
    g1 = download.download_genome("CP002779", test_temp_path)
    assert g.stat().st_mtime == g1.stat().st_mtime


def test_retries_download():
    a = [1]

    def call1(_: Path):
        if a[0] > 0:
            a[0] -= 1
            raise ConnectionRefusedError
        else:
            _.touch()

    download.retries_download(
        test_temp_path / "test_retries_download", call1, "", retry=3
    )
    a = 2
    download.retries_download(test_temp_path / ".." / "__init__.py", call1, "", retry=1)


@pytest.mark.xfail
def test_retries_download_fail():
    a = [1]

    def call1(_: Path):
        if a[0] > 0:
            a[0] -= 1
            raise ConnectionRefusedError
        else:
            _.touch()

    download.retries_download(
        test_temp_path / "test_retries_download_fail", call1, "", retry=1
    )
