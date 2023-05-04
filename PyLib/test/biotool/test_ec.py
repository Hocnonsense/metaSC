# -*- coding: utf-8 -*-
"""
 * @Date: 2023-05-04 11:21:46
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-04 21:59:10
 * @FilePath: /metaSC/PyLib/test/biotool/test_ec.py
 * @Description:
"""

from PyLib.biotool import ec
from PyLib.test import test_file_path


def test_check_title():
    ec_dat = test_file_path / "enzyme.dat"
    with ec_dat.open() as fi:
        is_title = [ec.check_title(i) for i in ec.raw_read_dat(fi)]
    assert is_title[0]
    assert sum(is_title[1:]) == 0


def test_check_entry_complete():
    ec_dat = test_file_path / "enzyme.dat"
    with ec_dat.open() as fi:
        is_entry_complete = [ec.check_entry_complete(i) for i in ec.raw_read_dat(fi)]
    assert not is_entry_complete[0]
    assert all(is_entry_complete[1:])


def test_parse_entry():
    ec_dat = test_file_path / "enzyme.dat"
    with ec_dat.open() as fi:
        for i, texts in zip(range(4), ec.raw_read_dat(fi)):
            pass
    entry_id, text_dict = ec.parse_entry_dict(texts)
    assert entry_id == "1.1.1.3"
    print()
    print(entry_id, text_dict)


def test_ec_entry():
    ec_dat = test_file_path / "enzyme.dat"
    with ec_dat.open() as fi:
        for texts in ec.raw_read_dat(fi):
            if ec.check_entry_complete(texts):
                ec_e = ec.EnzymeClassEntry.read_texts(texts)
                print()
                print(ec_e)


def test_ec_entry_str():
    ec_dat = test_file_path / "enzyme.dat"
    with ec_dat.open() as fi:
        for i, texts in zip(range(5), ec.raw_read_dat(fi)):
            pass
    ec_e = ec.EnzymeClassEntry.read_texts(texts)
    assert (
        str(ec_e)
        == "1.1.1.4 ((R,R)-butanediol dehydrogenase.) ['(R,R)-butane-2,3-diol + NAD(+) = (R)-acetoin + H(+) + NADH.']"
    )


def test_ec_entry_str_renamed():
    ec_dat = test_file_path / "enzyme.dat"
    with ec_dat.open() as fi:
        for i, texts in zip(range(6), ec.raw_read_dat(fi)):
            pass
    ec_e = ec.EnzymeClassEntry.read_texts(texts)
    assert str(ec_e) == "1.1.1.5 (Transferred entry) ['1.1.1.303', '1.1.1.304']"


def test_parse_id():
    assert ec.parse_id("1.1.1.3\n") == "1.1.1.3"


def test_bind_text():
    assert ec.bind_text("", "") == ""
    assert ec.bind_text("aa", "") == "aa"
    assert ec.bind_text("", "bb") == "bb"
    assert ec.bind_text("aa", "bb") == "aa bb"
    assert ec.bind_text("aa-", "bb") == "aa-bb"


def test_parse_de():
    assert (
        ec.parse_sentence(["Alcohol dehydrogenase.\n"])[0] == "Alcohol dehydrogenase."
    )
    assert (
        ec.parse_sentence(
            [
                "UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-\n",
                "alanyl-D-alanyl ligase.\n",
            ]
        )[0]
        == "UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanyl ligase."
    )
    assert (
        ec.parse_sentence(
            [
                "Transferred entry: 3.4.21.62, 3.4.21.63, 3.4.21.64, 3.4.21.65, 3.4.21.66\n",
                "and 3.4.21.67.\n",
            ]
        )[0]
        == "Transferred entry: 3.4.21.62, 3.4.21.63, 3.4.21.64, 3.4.21.65, 3.4.21.66 and 3.4.21.67."
    )


def test_parse_an():
    assert ec.parse_sentence(
        [
            "Terminal addition enzyme.\n",
            "Terminal deoxynucleotidyltransferase.\n",
            "Terminal deoxyribonucleotidyltransferase.\n",
            "Terminal transferase.\n",
        ]
    ) == [
        "Terminal addition enzyme.",
        "Terminal deoxynucleotidyltransferase.",
        "Terminal deoxyribonucleotidyltransferase.",
        "Terminal transferase.",
    ]


def test_parse_ca():
    assert ec.parse_sentence(
        [
            "L-malate + NAD(+) = oxaloacetate + NADH.\n",
            "2 ATP + NH(3) + CO(2) + H(2)O = 2 ADP + phosphate + carbamoyl\n",
            "phosphate.\n",
            "Cyclizes part of a 1,4-alpha-D-glucan chain by formation of a\n",
            "1,4-alpha-D-glucosidic bond.\n",
            "Cleavage of Leu-|-Xaa bond in angiotensinogen to generate\n",
            "angiotensin I.\n",
            "H2O + 2 NAD(+) + UDP-alpha-D-glucose = 3 H(+) + 2 NADH + UDP-alpha-D-\n",
            "glucuronate.\n",
        ]
    ) == [
        "L-malate + NAD(+) = oxaloacetate + NADH.",
        "2 ATP + NH(3) + CO(2) + H(2)O = 2 ADP + phosphate + carbamoyl phosphate.",
        "Cyclizes part of a 1,4-alpha-D-glucan chain by formation of a 1,4-alpha-D-glucosidic bond.",
        "Cleavage of Leu-|-Xaa bond in angiotensinogen to generate angiotensin I.",
        "H2O + 2 NAD(+) + UDP-alpha-D-glucose = 3 H(+) + 2 NADH + UDP-alpha-D-glucuronate.",
    ]


def test_parse_cc():
    assert ec.parse_comments(
        [
            "-!- The product spontaneously isomerizes to L-ascorbate.\n",
            "-!- Some members of this group oxidize only primary alcohols;\n",
            "    others act also on secondary alcohols.\n",
        ]
    ) == [
        "The product spontaneously isomerizes to L-ascorbate.",
        "Some members of this group oxidize only primary alcohols; others act also on secondary alcohols.",
    ]


def test_parse_dr():
    assert ec.parse_db_reference(
        [
            "P35497, DHSO1_YEAST;  Q07786, DHSO2_YEAST;  Q9Z9U1, DHSO_BACHD ;\n",
            "Q06004, DHSO_BACSU ;  Q02912, DHSO_BOMMO ;  Q00796, DHSO_HUMAN ;\n",
        ]
    ) == {
        "P35497": "DHSO1_YEAST",
        "Q07786": "DHSO2_YEAST",
        "Q9Z9U1": "DHSO_BACHD",
        "Q06004": "DHSO_BACSU",
        "Q02912": "DHSO_BOMMO",
        "Q00796": "DHSO_HUMAN",
    }


def test_extract_transferred_entries():
    assert ec.extract_transferred_entries("1.1.1.198, 1.1.1.227 and 1.1.1.228") == [
        "1.1.1.198",
        "1.1.1.227",
        "1.1.1.228",
    ]
