# -*- coding: utf-8 -*-
"""
 * @Date: 2022-11-12 20:31:16
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-04 10:38:09
 * @FilePath: /metaSC/PyLib/biotool/metabolic_overlap.py
 * @Description:
    extract from https://github.com/ericHester/metabolicOverlap/tree/master/scripts2

    for fasta in *.fna; do
        fna_to_rxn <$fasta | rxn_expandinfo >${fasta%.fna}.rxn
    done
    rxn_to_connections *.rxn | tee connections | lists_to_matrix -c3 >matrix

"""

import re
from pathlib import Path
from typing import Generator, NamedTuple, Literal, Union, overload
import pandas as pd


ModelSEEDDatabase = None


def set_seed_db(ModelSEEDDatabase_: Path = None):
    global ModelSEEDDatabase
    if ModelSEEDDatabase_:
        ModelSEEDDatabase = ModelSEEDDatabase_


def load_default_seed_db(ModelSEEDDatabase_: Path = None):
    return ModelSEEDDatabase_ or ModelSEEDDatabase


class Reaction(NamedTuple):
    id: str
    ins: list[str]
    outs: list[str]
    direction: Literal["=", ">", "<", "?"]

    def __hash__(self) -> int:
        return self.id.__hash__()

    @staticmethod
    def cleancompound(raw: str):
        """
        >>> cleancompound("cpd37299[0]")
        'cpd37299'
        """
        COMPOUND_RE = re.compile(r".*(cpd[0-9]+\[[^\]]*\])")
        if (match := re.match(COMPOUND_RE, raw)) is not None:
            return match.group(1)

    @classmethod
    def cleanequation(cls, id: str, raw: str, direction: Literal["=", ">", "<", "?"]):
        ins, outs = [
            [y for _y in x.split("+") if (y := cls.cleancompound(_y))]
            for x in re.split("<=>|=>|<=", raw)
        ]
        if direction == "<":  # Make sure it's always input > output
            direction = ">"
            ins, outs = outs, ins
        return cls(id, sorted(ins), sorted(outs), direction)


def load_clean_SEED_reactions(ModelSEEDDatabase):
    idx_id, idx_eq, idx_dir = 0, 6, 8
    filename = (
        load_default_seed_db(ModelSEEDDatabase) / "Biochemistry" / "reactions.tsv"
    )
    with open(filename) as file:
        next(file)
        return {
            row[idx_id]: Reaction.cleanequation(row[idx_id], row[idx_eq], row[idx_dir])
            for row in (line.rstrip("\n").split("\t") for line in file)
        }


def load_KO_modelseed_translations(ModelSEEDDatabase):
    KO_modelseed_translations_file = (
        load_default_seed_db(ModelSEEDDatabase)
        / "Biochemistry"
        / "Aliases"
        / "Provenance"
        / "Archived_Files"
        / "KO_modelseed_translations.csv"
    )
    ko_modelseed_translations: dict[str, list[str]] = {}
    with open(KO_modelseed_translations_file) as fi:
        for line in fi:
            values = line.strip().split()
            ko_modelseed_translations[values[0]] = values[1:]
    return ko_modelseed_translations


@overload
def genome_ko2rxn(
    kos: Generator[str, None, None],
    *,
    ko_modelseed_translations: dict[str, list[str]],
    clean_reactions: dict[str, Reaction],
) -> pd.Series:
    ...


@overload
def genome_ko2rxn(
    kos: Generator[str, None, None],
    *,
    ModelSEEDDatabase: Path,
) -> pd.Series:
    ...


def genome_ko2rxn(
    kos: Generator[str, None, None],
    ko_modelseed_translations: dict[str, list[str]] = None,
    clean_reactions: dict[str, Reaction] = None,
    ModelSEEDDatabase: Path = None,
):
    ko_modelseed_translations = (
        ko_modelseed_translations
        or load_KO_modelseed_translations(load_default_seed_db(ModelSEEDDatabase))
    )
    clean_reactions = clean_reactions or load_clean_SEED_reactions(
        load_default_seed_db(ModelSEEDDatabase)
    )

    return (
        pd.Series(kos)
        .apply(lambda i: ko_modelseed_translations.get(i, []))
        .explode()
        .sort_values()
        .dropna()
        .drop_duplicates()
        .apply(lambda i: clean_reactions[i])
    )


def genome_rxn2cpd(genome_rxn: Union[pd.Series, dict[str, list[Reaction]]]):
    return (
        pd.Series(genome_rxn)
        .explode()
        .rename("CPD")
        .pipe(
            lambda s: pd.concat(
                [
                    pd.DataFrame(
                        s.apply(
                            lambda rxn: {*rxn.ins, *rxn.outs}
                            if rxn.direction == "="
                            else set(rxn.ins)
                        ).explode()
                    ).assign(Substrate="Ins"),
                    pd.DataFrame(
                        s.apply(
                            lambda rxn: {*rxn.ins, *rxn.outs}
                            if rxn.direction == "="
                            else set(rxn.outs)
                        ).explode()
                    ).assign(Substrate="Outs"),
                ]
            )
            .reset_index()
            .rename({"index": "Genome"}, axis=1)
            .drop_duplicates()
        )
    )


def genome_rxn2mo(grxn):
    gcpd = genome_rxn2cpd(grxn)

    mo = pd.DataFrame(
        0,
        index=gcpd.groupby("Substrate")
        .get_group("Outs")["Genome"]
        .drop_duplicates()
        .rename("OutGenome"),
        columns=gcpd.groupby("Substrate")
        .get_group("Ins")["Genome"]
        .drop_duplicates()
        .rename("InGenome"),
    )
    for cpd, gcpd_i in gcpd.groupby("CPD"):
        if len(gcpd_i["Substrate"].drop_duplicates()) == 2:
            mo.loc[
                gcpd_i.groupby("Substrate").get_group("Outs")["Genome"],
                gcpd_i.groupby("Substrate").get_group("Ins")["Genome"],
            ] += 1

    return mo
