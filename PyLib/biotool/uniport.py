# -*- coding: utf-8 -*-
"""
 * @Date: 2022-09-23 23:07:44
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-24 12:12:48
 * @FilePath: /2021_09-MT10kSW/workflow/utils/PyLib/biotool/uniport.py
 * @Description:
from Bio.UniProt import GOA
"""

import gzip
import re
from pathlib import Path
from typing import TextIO, Type
from difflib import SequenceMatcher

import pandas as pd


def assign_unipprt_fa_features(
    headline=">ref_narH|tr|X5L138|X5L138_9MYCO Nitrate reductase subunit beta OS=Mycolicibacterium mageritense DSM 44476 = CIP 104973 OX=1209984 GN=narH PE=4 SV=1",
):
    features: dict[str, str] = {}
    if headline.startswith(">"):
        headline = headline[1:]
    headline = headline.strip()
    seqid, remain = headline.split(maxsplit=1)
    features["id"] = seqid
    key = "desc"
    while "=" in remain:
        value_, remain = remain.split("=", 1)
        while value_[-1] == " " or remain[0] == " ":
            remain_, remain = remain.split("=", 1)
            value_ = value_ + "=" + remain_
        features[key], key = value_.rsplit(maxsplit=1)
    features[key] = remain
    return features


class _UniProtTaxonomyDB:
    def __init__(self, taxonomy_file: Path):
        self.taxonomy_file = Path(taxonomy_file).expanduser().absolute()
        self.version = taxonomy_file.name
        assert self.taxonomy_file.is_file(), "please verify taxonomy_file"

    def load(self):
        try:
            with open(self.taxonomy_file) as ti:
                self._load_taxonomy(ti)
        except UnicodeDecodeError:
            with gzip.open(self.taxonomy_file) as ti:
                if self.version.endswith(".gz"):
                    self.version = self.version.rsplit(".gz", 1)[0]
                self._load_taxonomy(ti)

    def _load_taxonomy(self, ti: TextIO):
        self.taxonomy_db = pd.read_csv(ti, sep="\t").set_index("Taxon Id")
        if self.version.endswith(".tsv"):
            self.version = self.version.rsplit(".tsv", 1)[0]


class UniProtTaxonomy:
    _DB = _UniProtTaxonomyDB(Path(__file__))
    _cache: dict[int, "UniProtTaxonomy"] = {}

    def __init__(self, taxon_id: int):
        self.s = self._DB.taxonomy_db.loc[taxon_id, :]
        self.taxon_id = int(taxon_id)
        for k, v in self.s.to_dict().items():
            key = str(re.sub("[^a-zA-Z0-9_]", "_", k)).lower()
            # special optimization for [property:parent]
            if key not in self.__dir__():
                setattr(self, key, v)
        # print(f"adding {self} to database")
        self._cache[self.taxon_id] = self
        self._kids: list[UniProtTaxonomy] = None  # type: ignore

    def __new__(cls, taxon_id: int):
        if taxon_id not in cls._DB.taxonomy_db.index:
            return None
        if taxon_id in cls._cache:
            return cls._cache[taxon_id]
        return super().__new__(cls)

    def __repr__(self):
        rank = self.s["Rank"]
        name = self.s["Scientific name"]
        return f"[{rank}] {name}"

    def __str__(self):
        uid = self.taxon_id
        if self.is_root():
            lineage_info = "as root"
        else:
            cis_lineage = ", ".join(reversed(self.s["Lineage"].split(", ")))
            lineage_info = f"from [{cis_lineage}]"
        return f"{repr(self)} ({uid}) {lineage_info}"

    def is_root(self):
        return self.parent is None

    def get_parent(self, raw=False):
        parent = self.s["Parent"]
        if raw:
            return parent
        if parent != parent:
            print("found the root of this taxonomy database")
            return
        return type(self)(int(parent))

    parent = property(get_parent)

    @property
    def parents(self):
        parents: list[Type(self)] = []
        parent = self.parent
        while parent is not None and not parent.is_root():
            parents = [parent] + parents
            parent = parent.parent
        return parents

    @property
    def kids(self):
        if self._kids is None:
            try:
                kids = self._DB.taxonomy_db.groupby("Parent").get_group(self.taxon_id)
                self._kids = [type(self)(int(i)) for i in kids.index]
            except KeyError:
                self._kids = []
        return self._kids

    @classmethod
    def _init_wrapper(cls, up_tax_db, name="T"):
        return type(name, (cls,), dict(_DB=up_tax_db, _cache={}))

    @classmethod
    def init(cls, taxonomy_file: Path, name=None):
        uptaxon = _UniProtTaxonomyDB(taxonomy_file)
        uptaxon.load()
        return cls._init_wrapper(uptaxon, name or uptaxon.version)

    @classmethod
    def find(cls, somename, topn=10, raw=False):
        sth = SequenceMatcher()
        sth.set_seq2(somename)

        def match1(seq1):
            sth.set_seq1(seq1)
            sth.real_quick_ratio()
            return sth.ratio()

        top_hits = (
            cls._DB.taxonomy_db["Scientific name"]
            .apply(match1)
            .sort_values(ascending=False)
            .head(topn)
            .index
        )
        if raw:
            return UPTaxon._DB.taxonomy_db.loc[top_hits]
        else:
            return [cls(int(i)) for i in top_hits]


if __name__ == "__main__":
    UPTaxon: Type[UniProtTaxonomy] = UniProtTaxonomy.init(
        Path("~/Data/Database2")
        / "uniport/uniprot-taxonomy-2022.09.22-02.02.11.70.tsv.gz"
    )
    # UPTaxon._DB.taxonomy_db
    a = UPTaxon(555)
    a.parents
    # a.kids
    # UPTaxon._cache
    # type(a)
    UPTaxon.find("Nitroso")
