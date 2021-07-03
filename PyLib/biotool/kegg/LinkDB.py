# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-06 10:19:04
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-07-02 20:11:22
 * @FilePath: /2021_05-MT10kSW/Scripts/PyLib/biotool/kegg/LinkDB.py
 * @Description:
"""


from io import FileIO
from .query import load_brite, load_KEGG_module_raw
from .kmodule import KModule
from typing import List, Tuple, Dict, Union


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
