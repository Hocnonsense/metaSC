# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-06 10:19:04
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-06-06 14:26:33
 * @FilePath: /metaSC/PyLib/biotool/kegg/LinkDB.py
 * @Description:
"""


from io import FileIO
from typing import List, Tuple


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
