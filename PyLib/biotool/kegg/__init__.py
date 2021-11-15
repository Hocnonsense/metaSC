# -*- coding: utf-8 -*-
"""
 * @Date: 2021-03-28 12:36:23
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-11-15 20:57:26
 * @FilePath: /metaSC/PyLib/biotool/kegg/__init__.py
 * @Description:
"""
# noqa: F401
from .query import load_KEGG_module_raw, read_brite_json, load_brite, cached
from .kmodule import KModule, init_module
from .amino_acid_metabolism import AAm
from .LinkDB import path2tsv_iter, module_from_brite, load_KEGG_module, map_KO_dict, map_KO_substr
