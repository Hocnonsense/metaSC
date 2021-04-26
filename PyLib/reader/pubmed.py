# -*- coding: utf-8 -*-
"""
 * @Date: 2021-04-26 12:03:40
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-04-26 15:28:42
 * @FilePath: /metaSC/PyLib/reader/pubmed.py
 * @Description:
"""

import os
from io import FileIO, StringIO
from typing import Tuple
from PyLib.reader import path


def wrapStrFile(sth: FileIO, **kwargs):
    if isinstance(sth, str):
        if os.path.isfile(sth):
            sth = open(sth, **kwargs)
        else:
            sth = StringIO(sth)
    return sth


def add_pubmed_item(item: dict, k0: str, v: str):
    k1, v1 = '', v
    if v[-1] == ')':
        v, k1 = v[:-1].split('(')
        v1 = v.strip()
    elif v[-1] == ']':
        v, k1 = v[:-1].split('[')
        v1 = v.strip()
    item.setdefault(k0, {})[k1] = v1
    item.setdefault(k0, {})[''] = v


def pmiditer(text: FileIO = ''):
    text = wrapStrFile(text)
    while True:
        item = {}
        for line in text:
            line = line.rstrip('\n')
            if line[:6] == 'PMID- ':
                item['PMID'] = line[7:]
                break
        else:
            break
        for line in text:
            line = line.rstrip('\n')
            if not line:
                break
            if line[4:6] == '- ':
                k0, v = line[:4].strip(), ''
            v += line[6:]
            if not line.endswith(' '):
                add_pubmed_item(item, k0, v)
        yield item


def report(file, *keys: Tuple):
    for item in pmiditer(file):
        report_item = []
        for k in keys:
            if isinstance(k, str):
                k = k, ''
            k0, k1 = k
            report_item.append(item.get(k0, {'': ''})[k1])
        yield report_item


#from Bio import Entrez
#Entrez.email = 'hwrn.aou@sjtu.edu.cn'
#handle = Entrez.esummary(db='pubmed', id='33566237')
#record = Entrez.read(handle)
#record


if __name__ == '__main__':
    file = '/home/hwrn/Work/2021_02-MarianaMeta/Analyze/pubmed-27242695-set.nbib'
    dir = path.path(file)
    with open(dir.get('pubmed_set.2.tsv'), 'w') as fo:
        print('doi', 'Ti', 'JT', 'AB', 'DP', sep='\t', file=fo)
        for report_item in report(dir.get(), ('AID', 'doi'), 'TI', 'JT', 'AB', 'DP'):
            print('https://doi.org/', end='', file=fo)
            print(*report_item, sep='\t', file=fo)
