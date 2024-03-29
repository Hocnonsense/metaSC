# -*- coding: utf-8 -*-
"""
 * @Date: 2023-05-04 11:20:08
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-05 12:21:29
 * @FilePath: /metaSC/PyLib/biotool/ec.py
 * @Description: read ec text
"""
from enum import Enum
from typing import TextIO


def raw_read_dat(in_file: TextIO):
    texts: list[str] = []

    for line in in_file:
        if line.startswith("//"):
            yield texts
            texts = []
        else:
            texts.append(line)


def check_title(texts: list[str]):
    if {i[:5] for i in texts if i != "CC\n"} != {"CC   "}:
        return False
    if texts[0:5] != [
        "CC   -----------------------------------------------------------------------\n",
        "CC\n",
        "CC   ENZYME nomenclature database\n",
        "CC\n",
        "CC   -----------------------------------------------------------------------\n",
    ]:
        return False
    if (
        texts[-1]
        != "CC   -----------------------------------------------------------------------\n"
    ):
        return False
    return True


def check_entry_complete(texts: list[str]):
    ch = [i[:5] for i in texts]
    if ch[:2] != ["ID   ", "DE   "]:
        return False
    if "CA   " not in ch:
        if (
            len(set(ch)) == 2
            and EnzymeClassState.check_state(texts[1][5:]) != EnzymeClassState.valid
        ):
            return True
        return False
    for i in ch:
        if i not in ("ID   ", "DE   ", "AN   ", "CA   ", "CC   ", "PR   ", "DR   "):
            return False
    return True


def parse_entry_dict(texts: list[str]):
    entry_id = parse_id(texts[0][5:])
    text_dict: dict[str, list[str]] = {
        i: [] for i in ("DE", "AN", "CA", "CC", "PR", "DR")
    }
    for (k, v) in ((i[:2], i[5:]) for i in texts[1:]):
        text_dict.setdefault(k, []).append(v)
    return entry_id, text_dict


def parse_id(clean_text: str):
    return clean_text.rstrip()


def bind_text(text1: str, text2: str):
    if not text2:
        return text1
    if not text1:
        return text2
    if text1.endswith("-"):  # compounds
        return text1 + text2
    return text1 + " " + text2


def parse_sentence(clean_texts: list[str]):
    sentences: list[str] = []
    sentence = ""
    for i in clean_texts:
        sentence = bind_text(sentence, i.rstrip("\n"))
        if sentence.endswith("."):
            sentences.append(sentence)
            sentence = ""
    if sentence != "":
        sentences.append(sentence)
    return [i[:-1] for i in sentences]


def parse_ca(clean_texts: list[str]):
    if len(clean_texts) == 0:
        return []
    sentences: list[str] = []
    _i = 1
    li = f"({_i}) "
    sentence = ""
    if clean_texts[0].startswith(li):
        sentence = clean_texts[0][len(li) :].rstrip("\n")
        _i += 1
        li = f"({_i}) "
        for i in clean_texts[1:]:
            if i.startswith(li):
                sentences.append(sentence)
                _i += 1
                li = f"({_i}) "
                sentence = i[len(li) :].rstrip("\n")
            else:
                sentence = bind_text(sentence, i.rstrip("\n"))
    else:
        for i in clean_texts:
            sentence = bind_text(sentence, i.rstrip("\n"))
        print(sentence)
    if sentence != "":
        sentences.append(sentence)
    return [i[:-1] for i in sentences]


def parse_comments(clean_texts: list[str]):
    comments: list[str] = []
    comment = ""
    for i in clean_texts:
        if i.startswith("-!- "):
            if comment:
                comments.append(comment)
            comment = ""
        comment = bind_text(comment, i[4:].rstrip("\n"))
    if comment != "":
        comments.append(comment)
    return comments


def parse_db_reference(clean_texts: list[str]):
    return {
        k.strip(): v.strip()
        for k, v in (
            j.split(",")
            for i in (i.strip().split(";") for i in clean_texts)
            for j in i
            if j
        )
    }


class EnzymeClassState(Enum):
    valid = ""
    deleted = "Deleted entry"
    renumbered = "Transferred entry: "

    @classmethod
    def check_state(cls, DE: str):
        if DE.startswith(cls.deleted.value):
            return cls.deleted
        if DE.startswith(cls.renumbered.value):
            return cls.renumbered
        return cls.valid


def extract_transferred_entries(entries: str):
    return [i.strip(",") for i in entries.split() if i != "and"]


class EnzymeClassEntry:
    def __init__(
        self,
        ID: str,
        DE: str,
        AN: list[str] = None,
        CA: list[str] = None,
        CC: list[str] = None,
        PR: dict[str, str] = None,
        DR: dict[str, str] = None,
    ) -> None:
        self.ID = ID
        self.DE = DE
        self.AN = AN if AN is not None else []
        self.CA = CA if CA is not None else []
        self.CC = CC if CC is not None else []
        self.PR = PR if PR is not None else {}
        self.DR = DR if DR is not None else {}
        self.state = EnzymeClassState.check_state(self.DE)
        self.transfered = None
        if self.state == EnzymeClassState.renumbered:
            self.transfered = extract_transferred_entries(
                self.DE[len(self.state.value) :]
            )

    def __str__(self) -> str:
        if self.state == EnzymeClassState.deleted:
            return f"{self.ID} ({self.DE})"
        if self.state == EnzymeClassState.renumbered:
            return f"{self.ID} (Transferred entry) {self.transfered}"
        return f"{self.ID} ({self.DE}) {self.CA}"

    @classmethod
    def read_texts(cls, texts: list[str]) -> "EnzymeClassEntry":
        assert check_entry_complete(texts)
        entry_id, text_dict = parse_entry_dict(texts)
        return cls(
            ID=entry_id,
            DE=parse_sentence(text_dict["DE"])[0],
            AN=parse_sentence(text_dict["AN"]),
            CA=parse_ca(text_dict["CA"]),
            CC=parse_comments(text_dict["CC"]),
            PR=parse_db_reference(text_dict["PR"]),
            DR=parse_db_reference(text_dict["DR"]),
        )


class EnzymeClassDatabase:
    def __init__(self, dat: str, strict=True) -> None:
        self.dat = dat
        self.title: list[str] = []
        self.ecs = self.read_dat(strict=strict)

    def read_dat(self, strict=True):
        ecs: dict[str, EnzymeClassEntry] = {}
        with open(self.dat) as fi:
            text_iter = raw_read_dat(fi)
            if strict:
                title = next(text_iter)
                assert check_title(title)
                self.title = title
            for texts in text_iter:
                if check_entry_complete(texts):
                    e = EnzymeClassEntry.read_texts(texts)
                    ecs[e.ID] = e
                elif strict:
                    print(texts)
                    assert False
        if strict:
            for e in ecs.values():
                if e.state == EnzymeClassState.renumbered:
                    for i in e.transfered:
                        assert i in ecs
        return ecs

    def __contains__(self, key: str):
        return key in self.ecs

    def __getitem__(self, key: str) -> list[EnzymeClassEntry]:
        ec = self.ecs[key]
        if ec.state == EnzymeClassState.deleted:
            return []
        if ec.state == EnzymeClassState.renumbered:
            return [self.ecs.get(key, key) for key in ec.transfered]
        return [ec]
