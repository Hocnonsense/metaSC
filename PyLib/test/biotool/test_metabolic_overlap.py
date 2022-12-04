# -*- coding: utf-8 -*-
"""
 * @Date: 2022-12-04 09:21:36
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2022-12-04 09:32:36
 * @FilePath: /metaSC/PyLib/test/biotool/test_metabolic_overlap.py
 * @Description:
"""
# ###

from PyLib.biotool.metabolic_overlap import (
    Reaction,
    genome_ko2rxn,
    genome_rxn2mo,
    load_clean_SEED_reactions,
    load_default_seed_db,
    load_KO_modelseed_translations,
    set_seed_db,
)


def test_set_load_seed_db():
    set_seed_db("ModelSEEDDatabase")
    assert load_default_seed_db() == "ModelSEEDDatabase"
    assert load_default_seed_db("../ModelSEEDDatabase") != "ModelSEEDDatabase"
    set_seed_db("../ModelSEEDDatabase")
    assert load_default_seed_db() == "../ModelSEEDDatabase"


def test_Reaction():
    r1 = Reaction(
        id="rxn00012",
        ins=["cpd00076[0]"],
        outs=["cpd00027[0]", "cpd02298[0]"],
        direction="=",
    )
    r2 = Reaction.cleanequation(
        "rxn00012", "(2) cpd00076[0] <=> (1) cpd00027[0] + (1) cpd02298[0]", "="
    )
    assert r1 == r2
    assert r1.id == r2.id
    assert r1.ins == r2.ins
    assert r1.outs == r2.outs
    assert r1.direction == r2.direction


def demo():
    import pandas as pd

    ModelSEEDDatabase = "../ModelSEEDDatabase"
    ko_modelseed_translations = load_KO_modelseed_translations(ModelSEEDDatabase)
    clean_reactions = load_clean_SEED_reactions(ModelSEEDDatabase)

    genomeko = pd.read_csv("genomeko.csv", index_col=0)

    grxn: pd.Series = genomeko.apply(
        lambda kos: list(
            genome_ko2rxn(
                kos.pipe(lambda s: s.index[s > 0]),
                ko_modelseed_translations,
                clean_reactions,
            )
        )
    )
    mo = genome_rxn2mo(grxn)

    return mo
