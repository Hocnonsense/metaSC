from PyLib.PyLibTool.tmpPkl import Path as Path, TmpPkl as TmpPkl
from _typeshed import Incomplete
from typing import Dict, Type

def split_name(taxon__name: str): ...

class Taxon:
    UIDs: Dict
    names: Dict[str, Type['Taxon']]
    tax_n: Incomplete
    def __new__(cls, name): ...
    taxon: Incomplete
    name: Incomplete
    parent: Incomplete
    kids: Incomplete
    def __init__(self, name) -> None: ...
    def get_UID(self): ...
    def set_UID(self, UID: str): ...
    UID: Incomplete
    def get_super(self, taxon): ...
    @classmethod
    def add(cls, long_taxonomy: str): ...
    @classmethod
    def attach_UID(cls, long_taxonomy: str, UID): ...

def read_taxonomy(CHECKM_DIR) -> None: ...
def read_metadata(CHECKM_DIR) -> None: ...

CHECKM_DIR: Incomplete
