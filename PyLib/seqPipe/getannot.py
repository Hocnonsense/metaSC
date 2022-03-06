# -*- coding: utf-8 -*-
"""
 * @Date: 2022-03-01 10:55:42
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-03-04 15:16:15
 * @FilePath: /metaSC/PyLib/seqPipe/getannot.py
 * @Description:
"""

from collections import OrderedDict
import os
import re
from typing import Dict, Iterable, List, Set, TextIO, Tuple, Union

from PyLib.PyLibTool.file_info import verbose_import
from PyLib.reader.iters import emapper_iter, read_table

logger = verbose_import(__name__, __doc__)


def drop_gene(gene: str, subsets: Tuple[Union[Set, List, Dict], bool]) -> bool:
    subset, subset_is_contig = subsets
    if subset:
        if subset_is_contig:
            return gene.rsplit("_", 1)[0] not in subset
        else:
            return gene not in subset
    return False


def GhostKOALA_iter(text: TextIO) -> Iterable[Tuple[str, str]]:
    for values in read_table(text):
        if len(values) > 1:
            gene, ko = values[0:2]
            yield gene, ko


def KofamKOALA_iter(text: TextIO) -> Iterable[Tuple[str, str]]:
    for values in read_table(text):
        if len(values) > 1:
            gene, ko = values[1:3]
            yield gene, ko


def eggnog_iter(text: TextIO) -> Iterable[Tuple[str, str]]:
    """Only report the first match:
    >>> yield gene, kos.split(",")[0][3:]"""
    i_KEGG_ko = 11
    for values in emapper_iter(text):
        kos = values[i_KEGG_ko]
        if kos:
            gene = values[0]
            yield gene, kos.split(",")[0][3:]


## collect gene KO
formats_func = OrderedDict(
    ghost=GhostKOALA_iter,
    kofam=KofamKOALA_iter,
    eggnog=eggnog_iter,
)


def get_gene_KOs(
    ann_files: List[str], subsets: Tuple[Union[Set, List, Dict], bool] = None
) -> Dict[str, str]:
    """Only keep the first match:
    >>> gene_KOs.setdefault(gene, ko)"""
    gene_KOs: Dict[str, str] = {}
    for format, file in zip(formats_func, ann_files):
        if not os.path.exists(file):
            logger.warning(f"{format} file does not exist, skip")
            continue
        logger.info(f"reading {file}")

        with open(file) as file_in:
            for gene, ko in formats_func[format](file_in):
                if subsets and drop_gene(gene, subsets):
                    continue
                gene_KOs.setdefault(gene, ko)
        logger.warning(f"{len(gene_KOs)} genes annotated...")

    return gene_KOs


def infer_annot_filename(pattern: str):
    ann_files = ["", "", ""]

    in_dir = ""
    index = pattern.rfind("/") + 1
    if index:
        in_dir, pattern = pattern[:index], pattern[index:]
        if not os.path.exists(in_dir):
            logger.warning("illigal pattern, ignore")
            pattern = ""
    # find avaiable annotation file
    pattern_re = re.compile(pattern)
    for file in sorted(os.listdir(in_dir)):
        if pattern_re.search(file):
            for i, source in enumerate(formats_func):
                if source in file.lower():
                    ann_files[i] = os.path.join(in_dir, file)

    if not any(ann_files):
        logger.fatal("illigal pattern, please check!")
        raise FileNotFoundError(f"pattren '{pattern}' donot match any file!")

    return ann_files
