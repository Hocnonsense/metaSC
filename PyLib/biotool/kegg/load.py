# -*- coding: utf-8 -*-
"""
 * @Date: 2023-02-08 11:23:52
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-02-08 11:28:09
 * @FilePath: /2022_09-M_mem/workflow/utils/libs/metaSC/PyLib/biotool/kegg/load.py
 * @Description:
"""
# """

import pandas as pd

from pathlib import Path
from .LinkDB import module_from_brite


def load_ko00002(KEGG_DIR="./KEGG_DB"):
    """Database may be download from KEGG, including the file of module and description (ko00002.json)"""
    import pandas as pd

    KEGG_DIR = Path(KEGG_DIR).expanduser()

    KEGG_DIR.mkdir(parents=True, exist_ok=True)

    module_levels, modules = module_from_brite(
        "br:ko00002",
        KEGG_DIR / "brite" / "ko00002.json",
        KEGG_DIR / "module",
    )
    module_levels_ = pd.DataFrame(
        module_levels, columns=["A", "B", "C", "entry", "name"]
    )
    module_levels_.index = module_levels_["entry"]  # type: ignore
    return module_levels_, dict(modules)


def load_entry(KEGG_DIR="./KEGG_DB"):
    module_levels, modules = load_ko00002(KEGG_DIR)
    entry2ko = pd.concat(
        [
            pd.DataFrame({"entry": entry, "KO": module.kos})
            for entry, module in modules.items()
        ]
    )
    # module_levels.to_csv("data/module_levels.tsv", sep="\t", index=False)
    # entry2ko.to_csv("data/entry2ko.tsv", sep="\t", index=False)
    return module_levels, entry2ko


def get_gmodule(genomeko: pd.DataFrame, KEGG_DIR="./KEGG_DB"):
    _, modules = load_ko00002(KEGG_DIR)
    gmodule_ = (
        genomeko.apply(lambda x: x[x > 0].index, axis=0)
        .apply(
            lambda x: {
                mname: module.completeness(x) for mname, module in modules.items()
            }
        )
        .apply(lambda x: pd.Series(x))
    )
    gmodule = gmodule_.T[gmodule_.apply(sum, 0) > 0]

    return gmodule


def demo2():
    module_levels, modules = load_ko00002()
    ko_abd = {
        "K19746": 0.0,
        "K19744": 1.0,
        "K12658": 2.0,
        "K21060": 3.0,
        "K22549": 4.0,
        "K21061": 5.0,
        "K22550": 6.0,
        "K21062": 7.0,
        "K13877": 8.0,
    }
    print(
        modules["M00947"].abundance(ko_abd),
        modules["M00947"].completeness(ko_abd),
        modules["M00947"].completeness(
            {ko: abd for ko, abd in ko_abd.items() if abd > 0}
        ),
    )
    print(
        modules["M00948"].abundance(ko_abd),
        modules["M00948"].completeness(ko_abd),
    )