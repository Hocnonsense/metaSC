# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-14 18:41:24
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-02-08 11:26:21
 * @FilePath: /2022_09-M_mem/workflow/utils/libs/metaSC/PyLib/biotool/kegg/query.py
 * @Description:
"""

import os
import json
from typing import Callable, List, Dict, TextIO, Union

from Bio.KEGG import REST


def load_KEGG_module_raw(
    source: Union[str, TextIO], cache_path: str = ""
) -> Dict[str, List[str]]:
    if isinstance(source, str):
        if cache_path:
            source_ = cached(_load_module_txt, source, cache_path)
        else:
            source_ = _load_module_txt(source)
    with source_ as file:  # type: ignore
        raw_module: Dict[str, List[str]] = {}
        listv = []
        for line in file:
            if line.startswith("///"):
                break
            is_k, v = line[:12].strip(), line[12:].strip()
            if is_k:
                k = is_k
                listv = raw_module.setdefault(k, [])
            if v:
                listv.append(v)

    return raw_module


def read_brite_json(brite_doc):
    name, children = brite_doc["name"], brite_doc.get("children", "")
    if not children:
        return name.split(maxsplit=1)  # type: ignore
    return name, dict((read_brite_json(children_doc) for children_doc in children))


def load_brite(source: Union[str, TextIO], cache_path: str = ""):
    if isinstance(source, str):
        if cache_path:
            source_ = cached(_load_brite_json, source, cache_path)
        else:
            source_ = _load_brite_json(source)
    with source_ as json_in:  # type: ignore
        brite_doc: Dict[str, Union[str, List[Dict]]] = json.loads(json_in.read())
        return read_brite_json(brite_doc)


def _load_module_txt(source: str) -> TextIO:
    """Resolve source to file or download it directly."""
    if os.path.isfile(source):
        return open(source)
    if source.startswith("M"):
        return REST.kegg_get(source)
    raise FileNotFoundError("not a file nor a right module number")


def _load_brite_json(source: str) -> TextIO:
    """Resolve source to file or download it directly."""
    if os.path.isfile(source):
        return open(source)
    if source.startswith("br:ko"):
        return REST.kegg_get(source, "json")
    raise FileNotFoundError("not a file nor a right br number")


def cached(
    func: Callable[[str], TextIO], source: str, filename: Union[None, str] = None
):
    # case 1: source or filename is file:
    if os.path.isfile(source):
        return func(source)

    filename = os.path.abspath(os.path.expanduser(filename or "./"))
    if os.path.isfile(filename):
        return func(filename)
    # else source MUST be KEGG id.
    # case 2: filename is directory:
    if os.path.isdir(filename):
        for file in os.listdir(filename):
            if source in file:
                return func(os.path.join(filename, file))
        # case 3: no cached, record file from func and return read(file)
        filename = os.path.join(filename, source)
    # case 4: no cached, will write to filename
    with open(filename, "w") as file_out:
        file_in: TextIO = func(source)
        file_out.write(file_in.read())
    return open(os.path.join(filename))
