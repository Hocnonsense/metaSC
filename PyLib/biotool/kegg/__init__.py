# -*- coding: utf-8 -*-
"""
 * @Date: 2021-03-28 12:36:23
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-02-08 11:32:16
 * @FilePath: /2022_09-M_mem/workflow/utils/libs/metaSC/PyLib/biotool/kegg/__init__.py
 * @Description:
"""
from .query import (
    load_KEGG_module_raw,
    read_brite_json,
    load_brite,
    cached,
)
from .kmodule import KModule, init_module
from .amino_acid_metabolism import AAm
from .LinkDB import (
    path2tsv_iter,
    module_from_brite,
    load_KEGG_module,
    map_KO_dict,
    map_KO_substr,
)
from .load import load_ko00002, load_entry, get_gmodule
