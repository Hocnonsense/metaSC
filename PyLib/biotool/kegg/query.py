# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-14 18:41:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-07-02 16:59:23
 * @FilePath: /metaSC/PyLib/biotool/kegg/query.py
 * @Description:
"""

import os
from io import FileIO
import json
from typing import Callable, List, Dict, Tuple, Union

from Bio.KEGG import REST

from PyLib.biotool.kegg.kmodule import KModule


def load_KEGG_module_raw(source: Union[str, FileIO], cache_path: str = ''
                         ) -> Dict[str, List[str]]:
    if isinstance(source, str):
        if cache_path:
            source = cached(_load_module_txt, source, cache_path)
        else:
            source = _load_module_txt(source)
    with source as file:
        raw_module: Dict[str, List[str]] = {}
        listv = []
        for line in file:
            if line.startswith('///'):
                break
            is_k, v = line[:12].strip(), line[12:].strip()
            if is_k:
                k = is_k
                listv = raw_module.setdefault(k, [])
            listv.append(v)

    return raw_module


def read_brite_json(brite_doc: Dict[str, Union[str, List[Dict]]]):
    name, children = brite_doc['name'], brite_doc.get('children', '')
    if not children:
        return name.split(maxsplit=1)
    return name, dict((read_brite_json(children_doc) for children_doc in children))


def load_brite(source: Union[str, FileIO], cache_path: str = ''
               ) -> Tuple[str, Dict[str, Dict[str, Dict[str, Dict[str, str]]]]]:
    if isinstance(source, str):
        if cache_path:
            source = cached(_load_brite_json, source, cache_path)
        else:
            source = _load_brite_json(source)
    with source as json_in:
        brite_doc: Dict[str, List[Dict]] = json.loads(json_in.read())
        return read_brite_json(brite_doc)


def _load_module_txt(source: str) -> FileIO:
    """Resolve source to file or download it directly."""
    if os.path.isfile(source):
        return open(source)
    if source.startswith('M'):
        return REST.kegg_get(source)
    raise FileNotFoundError('not a file nor a right module number')


def _load_brite_json(source: str) -> FileIO:
    """Resolve source to file or download it directly."""
    if os.path.isfile(source):
        return open(source)
    if source.startswith('br:ko'):
        return REST.kegg_get(source, 'json')
    raise FileNotFoundError('not a file nor a right br number')


def cached(func: Callable[[str], FileIO], source: str,
           filename: Union[None, str] = None):
    # case 1: source or filename is file:
    if os.path.isfile(source):
        return func(source)

    filename = os.path.abspath(os.path.expanduser(filename or './'))
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
    with open(filename, 'w') as file_out:
        file_in: FileIO = func(source)
        file_out.write(file_in.read())
    return open(os.path.join(filename))


def demo1():
    KEGG_DIR = '~/Work/2021_05-MT10kSW/00_data/KEGG/ko00002'
    ko00002 = load_brite('br:ko00002', os.path.join(KEGG_DIR, 'ko00002.json'))[1]

    module_name = list(list(list(list(
        ko00002.values()
    )[0].values())[0].values())[0].keys())[0]
    raw_module = load_KEGG_module_raw(module_name, KEGG_DIR)
    entry: str = raw_module['ENTRY'][0].split()[0]
    module = KModule(''.join(raw_module['DEFINITION']),
                     additional_info=''.join(raw_module['NAME']))
    print(entry)
    return module
