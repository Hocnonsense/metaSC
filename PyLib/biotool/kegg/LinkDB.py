# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-06 10:19:04
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-11-15 20:55:11
 * @FilePath: /metaSC/PyLib/biotool/kegg/LinkDB.py
 * @Description:
"""

import re
from io import FileIO
from typing import Dict, List, Tuple, Union

from .kmodule import KModule
from .query import load_brite, load_KEGG_module_raw


def path2tsv_iter(text: FileIO) -> Tuple[str, List[str], str, List[str]]:
    """
     * @description: read word from
        https://www.genome.jp/dbget-bin/get_linkdb?-t+orthology+path:ko00910
        and change it to tsv
     * @param text: string starts with
        ID                   Definition
        ----------------------------------------------------------------------------------------------------
     * @return {*}
        ID, gene, desc, EC
    """
    for line in text:
        if line.startswith('K') and set(line[1:6]).issubset(set('1234567890')):
            ID, values = line.split(' ', 1)
            values = values.strip()
            genes, values = values.split('; ', 1)
            genes = genes.split(', ')
            values = values.strip()
            desc, EC = (values.split('[EC:', 1) + [''])[:2]
            desc: str = desc.strip()
            EC = EC.split(']', 1)[0].split()
            yield ID, genes, desc, EC


def module_from_brite(source: Union[str, FileIO], brite_path: str = '', module_path: str = ''
                      ) -> Tuple[List[str], List[Tuple[str, KModule]]]:
    _, brite = load_brite(source, brite_path)
    module_levels = []
    modules: List[Tuple[str, KModule]] = []
    for modules1_name, modules1 in brite.items():
        for metabolism_name, metabolism in modules1.items():
            for metabolism_name2, metabolism2 in metabolism.items():
                for entry, name in metabolism2.items():
                    module_levels.append((modules1_name, metabolism_name,
                                          metabolism_name2, entry, name))
                    raw_module = load_KEGG_module_raw(entry, module_path)
                    module = KModule(''.join(raw_module['DEFINITION']),
                                     additional_info=''.join(raw_module['NAME']))
                    modules.append((entry, module))
    return module_levels, modules


def load_KEGG_module(module_name, KEGG_DIR):
    raw_module = load_KEGG_module_raw(module_name, KEGG_DIR)
    return KModule(''.join(raw_module['DEFINITION']), additional_info=raw_module)


def map_KO_dict(KO_str: str, ko_match: Dict[str, float]) -> str:
    """ Map KO abundance to a KO str
    """
    KO_PATTERN = re.compile(r"K\d{5}")
    for match in re.finditer(KO_PATTERN, KO_str):
        KO = match.group()
        if KO in ko_match:
            KO_value = ko_match[KO]
            KO_value = "{value:.{len_value}f}".format(
                value=KO_value, len_value=5-len(str(KO_value // 1).split(".")[0]))
        else:
            KO_value = " NaN. "
        KO_str = KO_str.replace(KO, KO_value)
    return KO_str


def map_KO_substr(KO_str: str, KO_substr: str) -> str:
    """ Map a str of KO to another KO str
        The KO in the first string but not in the second string will be masked
    """
    KO_PATTERN = re.compile(r"K\d{5}")
    for match in re.finditer(KO_PATTERN, KO_str):
        KO = match.group()
        if KO in KO_substr:
            KO_value = KO
        else:
            KO_value = "      "
        KO_str = KO_str.replace(KO, KO_value)
    return KO_str
