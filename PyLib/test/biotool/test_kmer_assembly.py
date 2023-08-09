# -*- coding: utf-8 -*-
"""
 * @Date: 2023-08-08 20:18:58
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-09 15:09:11
 * @FilePath: /metaSC/PyLib/test/biotool/test_kmer_assembly.py
 * @Description:
"""

from PyLib.biotool.kmer_assembly import KmerAssembly
from PyLib.test import test_temp_path


def test_kmer_assembly():
    test_fa = (
        ">k7_1\nTTCTCGATAGGCTT\n"
        ">k7_2\nGGCTTCTATCCCGC\n"
        ">k7_3'\nAAGCCTAGGTCTCA\n"
        ">k7_4\nCTAGACCTATTCTC\n"
        ">k7_5c\nCAAGACCTACAAGA\n"
        ">k7_6\nGCTATTCGAGCAGC\n"
    )
    file = test_temp_path / "test_kmer_assembly.fa"
    with open(file, "w") as fo:
        fo.write(test_fa)

    ka = KmerAssembly(file, 5)
    ka.get_kmer()
    ka.filter_ends()

    assert set(ka.kmers) == {"AAGCC", "GAGAA", "CAAGA"}
    assert [ka.kmers[kmer].report_dict() for kmer in sorted(ka.kmers)] == [
        {
            "kmer": "AAGCC",
            "starts": "k7_1,k7_3'",
            "ends": "k7_2",
            "circle": "",
            "single": "k7_2,k7_3'",
        },
        {
            "kmer": "CAAGA",
            "starts": "",
            "ends": "",
            "circle": "k7_5c",
            "single": "",
        },
        {
            "kmer": "GAGAA",
            "starts": "k7_4",
            "ends": "k7_1",
            "circle": "",
            "single": "k7_4",
        },
    ]
