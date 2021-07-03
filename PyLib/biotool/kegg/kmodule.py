# -*- coding: utf-8 -*-
"""
 * @Date: 2020-07-01 00:29:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-06-15 12:11:38
 * @FilePath: /metaSC/PyLib/biotool/kegg/kegg.py
 * @Description:
"""
import pickle
from io import StringIO
from typing import Dict, Hashable, List


class KModule():
    """
        @description:
            KEGG module units
            It contains the functions in elements.

            if self.ko, then assert not self.elements, vice versa.
    """
    def __init__(self, express="", additional_info=""):
        # print(express)
        self.is_chain: bool = True
        self.elements: List[KModule] = []
        self.ko: str = ""
        self.additional_info: str = additional_info

        self.__calculate(express)

    def __calculate(self, express: str):
        """Calculate express"""
        if not express:
            return
        # print(express, "<")
        p = len(express) - 1
        if p == 5:  # ko
            self.ko = express
            return
        no_comma = []
        while p >= 0:
            c = express[p]
            if c in "0123456789":  # get a brief KO numbner
                no_comma.append(KModule(express[p - 5: p + 1]))
                p -= 5
            elif c in ")":  # you are trapped into a sub element
                last_p = p
                bracket_stack = 1
                plb = express.rfind("(", 0, p)
                prb = express.rfind(")", 0, p)
                while bracket_stack:
                    if plb < prb:
                        bracket_stack += 1
                        p = prb
                        prb = express.rfind(")", 0, p)
                    else:
                        bracket_stack -= 1
                        p = plb
                        plb = express.rfind("(", 0, p)
                no_comma.append(KModule(express[p + 1:last_p]))
            elif c in ",":
                if self.is_chain:
                    self.is_chain = False
                # Stupid KEGG think "(K09880,K08965 K08966)" is equal to "(K09880,(K08965 K08966))".
                no_comma.reverse()  # "," is in this express. Reverse "no_comma"
                self.elements.append(KModule.from_list(no_comma))
                no_comma = []
            p -= 1
        if self.elements:  # "," is in this express. Reverse "no_comma"
            no_comma.reverse()
            self.elements.append(KModule.from_list(no_comma))
        else:  # no "," is in this express, so just add it.
            self.elements = no_comma
        # print(*self.elements)
        self.elements.reverse()

    def abundance(self, ko_match: Dict[str, float]):
        return sum(ko_match.get(ko, 0) for ko in self.list_ko())

    def completeness(self, ko_match: Hashable) -> float:
        """Complessness of given match, ko is its dict"""
        count = 0
        if self.ko:
            return 1 if self.ko in ko_match else 0
        # multipy elements
        if self.is_chain:
            for element in self.elements:
                count += element.completeness(ko_match)
            return count / len(self.elements)
        # self.is_chain is False
        return max([element.completeness(ko_match) for element in self.elements])

    def list_ko(self) -> List[str]:
        if self.ko:
            return [self.ko]
        return [ko for element in self.elements for ko in element.list_ko()]

    def __len__(self):
        return sum([len(e) for e in self.elements]) if self.elements else 1

    def __getitem__(self, key):
        """Return a way contains it"""
        if key == self.ko:
            return [self.ko], [0]
        for e in self.elements:
            e_key, i_key = e[key]
            if i_key != -1:
                if self.is_chain:  # all elements is important
                    e_key = [e_key if e_chain is e else e_chain for e_chain in self.elements]
                    i_key = [self.elements.index(e)] + i_key
                return e_key, i_key
        return [], -1

    def __str__(self):
        if self.ko:
            return self.ko
        sep = " " if self.is_chain else ","
        return sep.join([kid.ko if kid.ko else "(" + str(kid) + ")" for kid in self.elements])

    @classmethod
    def from_list(cls, no_comma, is_chain=True, additional_info=""):
        """New Element by a list"""
        if len(no_comma) == 1:
            e = no_comma[0]
        else:
            e = KModule()
            e.elements = no_comma
            e.is_chain = is_chain
            e.additional_info = additional_info
        return e


def init_module(module_str):
    """
        @description: init modules
        @param module_str:
            string or IO of modules from KEGG module
            e.g. Amino_acid_metabolism
        @return module{metabolism: {Entry: KModule(Definition)}}
    """
    if isinstance(module_str, str):
        module_str = StringIO(module_str)

    module = {}
    title = ""
    line = module_str.readline()
    while line:
        if line.strip()[0:5] == "Entry":
            Entry = line.split()[1]
            line = module_str.readline()
            Name = " ".join(line.strip().split()[1:])
            # Name = Name[: Name.find(",")]
            line = module_str.readline()
            express = " ".join(line.strip().split()[1:])
            module[title][Entry] = KModule(express, Name)
        else:  # title
            title = line.strip()[:-1]
            if title:
                module[title] = {}
        line = module_str.readline()

    return module


def _load_module(module_name):
    """
        @description: load module by pickle files or text files
    """
    try:
        pin = open(module_name + ".pickle", "rb")
        module = pickle.load(pin)
    except (OSError, FileNotFoundError, EOFError) as exc:
        if exc.errno != 36:  # file too long :)
            raise
        module = init_module(module_name)
        # with open(module_name + ".pickle", "wb") as pout:
        #     pickle.dump(module, pout)
    return module
