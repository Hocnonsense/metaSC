# -*- coding: utf-8 -*-
"""
 * @Date: 2021-05-04 16:41:43
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-18 11:13:00
 * @FilePath: /metaSC/PyLib/biotool/checkm/UID_taxon.py
 * @Description:
    find and match UID and taxon.
"""

import os
from typing import Dict, Type, Union

from PyLib.PyLibTool.tmpPkl import TmpPkl


def split_name(taxon__name: str):
    return taxon__name.split('__')


class Taxon:
    UIDs = {}
    names: Dict[str, Type['Taxon']] = {}
    tax_n = ('k', 'p', 'c', 'o', 'f', 'g', 's')

    def __init__(self, name) -> None:
        if name in Taxon.names:
            self = Taxon.names[name]
            return
        Taxon.names[name] = self
        self.taxon: str = name[0]
        self.name = name
        self.__UID = ''
        self.parent: Taxon = None
        self.kids = set()

    def __str__(self) -> str:
        return self.name

    def get_UID(self):
        return self.__UID

    def set_UID(self, UID: str):
        #assert not self.__UID  # init only once
        self.__UID = UID
        Taxon.UIDs[UID] = self
    UID = property(fget=get_UID, fset=set_UID)

    def get_super(self, taxon) -> Union[None, Type['Taxon']]:
        parent = self
        while parent and parent.taxon != taxon:
            parent = parent.parent
        return parent

    @classmethod
    def add(cls, long_taxonomy: str):
        Flag = False

        taxonomy_iter = reversed(long_taxonomy.split(';'))
        taxonomy: str = next(taxonomy_iter)
        while taxonomy.endswith('__'):
            taxonomy: str = next(taxonomy_iter)
        if taxonomy in Taxon.names:
            return

        kid = cls(taxonomy)
        for taxonomy in taxonomy_iter:
            if taxonomy in Taxon.names:
                parent = cls.names[taxonomy]
                Flag = True
            else:
                parent = cls(taxonomy)

            kid.parent = parent
            kid = parent
            parent.kids.add(kid)

            if Flag:
                break

    @classmethod
    def attach_UID(cls, long_taxonomy: str, UID: str):
        long_taxonomy: list = long_taxonomy.split(';')
        for taxonomy in long_taxonomy:
            if taxonomy not in Taxon.names:
                Taxon.add(taxonomy)
            Taxon.names[taxonomy].UID = UID


def read_taxonomy(CHECKM_DIR):
    #taxonomy: Dict[str, Dict[str, List[str]]] = {
    #    taxon: {} for taxon in ('k', 'p', 'c', 'o', 'f', 'g', 's')}
    with open(os.path.join(CHECKM_DIR, 'genome_tree', 'genome_tree.taxonomy.tsv')) as ti:
        for line in ti:
            values = line.strip().split('\t')
            if values:
                _, taxons = values[0], values[1]
                Taxon.add(taxons)


def read_metadata(CHECKM_DIR):
    #UID_taxon: Dict[str, List[str]] = {}
    with open(os.path.join(CHECKM_DIR, 'genome_tree', 'genome_tree.metadata.tsv')) as mi:
        header = mi.readline()
        assert header.startswith('UID')
        for line in mi:
            values = line.strip().split('\t')
            if values and values[2]:
                UID, taxons = values[0], values[2]  #.split(';')[-1]
                Taxon.attach_UID(taxons, UID)


with TmpPkl('UID_2015_01_16.pickle', desc=__doc__,
            force_rewrite=True,
            situ=__file__) as tmp:
    if tmp.force_rewrite:
        CHECKM_DIR = os.path.expanduser('~/.checkm')
        read_taxonomy(CHECKM_DIR)
        read_metadata(CHECKM_DIR)
        Taxon.add('root')
        tmp.last_results = Taxon.names, Taxon.UIDs

    Taxon.names, Taxon.UIDs = tmp.last_results
