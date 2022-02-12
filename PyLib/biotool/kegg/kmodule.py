# -*- coding: utf-8 -*-
"""
 * @Date: 2020-07-01 00:29:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-09 10:26:07
 * @FilePath: /metaSC/PyLib/biotool/kegg/kmodule.py
 * @Description:
"""
import pickle
from io import StringIO
from typing import Dict, List, Set, Union


class KModule:
    """
    @description:
        KEGG module units
        It contains the functions in elements.

        if self.ko, then assert not self.elements, vice versa.
    """

    def __init__(self, express="", additional_info=""):
        # print(express)
        self._is_chain: bool = True
        self._elements: List[KModule] = []
        self._ko: str = ""
        self.additional_info: str = additional_info

        self.__calculate(express)

    def __calculate(self, express: str):
        """Calculate express"""
        if not express:
            return
        # print(express, "<")
        p = len(express) - 1
        if p == 5:  # ko
            self._ko = express
            return
        no_comma = []
        while p >= 0:
            c = express[p]
            if c in "0123456789":  # get a brief KO numbner
                no_comma.append(KModule(express[p - 5 : p + 1]))
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
                no_comma.append(KModule(express[p + 1 : last_p]))
            elif c in ",":
                if self._is_chain:
                    self._is_chain = False
                # Stupid KEGG think "(K09880,K08965 K08966)" is equal to "(K09880,(K08965 K08966))".
                no_comma.reverse()  # "," is in this express. Reverse "no_comma"
                self._elements.append(KModule.from_list(no_comma))
                no_comma = []
            p -= 1
        if self._elements:  # "," is in this express. Reverse "no_comma"
            no_comma.reverse()
            self._elements.append(KModule.from_list(no_comma))
        else:  # no "," is in this express, so just add it.
            self._elements = no_comma
        # print(*self.elements)
        self._elements.reverse()

    def list_ko(self) -> List[str]:
        if self._ko:
            return [self._ko]
        return [ko for element in self._elements for ko in element.list_ko()]

    kos = property(fget=list_ko)

    def __len__(self):
        return sum([len(e) for e in self._elements]) if self._elements else 1

    def __getitem__(self, key):
        """Return a way contains it"""
        if key == self._ko:
            return [self._ko], [0]
        for e in self._elements:
            e_key, i_key = e[key]
            if i_key != -1:
                if self._is_chain:  # all elements is important
                    e_key = [
                        e_key if e_chain is e else e_chain for e_chain in self._elements
                    ]
                    i_key = [self._elements.index(e)] + i_key
                return e_key, i_key
        return [], -1

    def __str__(self):
        if self._ko:
            return self._ko
        sep = " " if self._is_chain else ","
        return sep.join(
            [kid._ko if kid._ko else "(" + str(kid) + ")" for kid in self._elements]
        )

    @classmethod
    def from_list(cls, no_comma, is_chain=True, additional_info=""):
        """New Element by a list"""
        if len(no_comma) == 1:
            e = no_comma[0]
        else:
            e = KModule()
            e._elements = no_comma
            e._is_chain = is_chain
            e.additional_info = additional_info
        return e

    def all_paths(self, ko_match: Union[List, Dict, Set] = None) -> List[str]:
        """
        module.all_paths():
            return all potential metabolism paths with KO
        module.all_paths(ko_match):
            return all available metabolism paths if KO be found in ko_match
        """
        if self._is_chain:
            if self._ko:
                if ko_match is None or self._ko in ko_match:
                    return [self._ko]
                return []  # this KO is not in list/dict/set, should not be detected
            last_paths = [""]
            for kid in self._elements:
                paths = kid.all_paths(ko_match)
                if not paths:
                    return []
                last_paths = [
                    f"{last_path} {path}" if last_path else path
                    for last_path in last_paths
                    for path in paths
                ]
            return last_paths
        else:
            paths = [path for kid in self._elements for path in kid.all_paths(ko_match)]
            if paths == []:
                return []
            len_path = max(len(path) for path in paths)
            paths = sorted(
                {
                    "[{path:^{len_path}}]".format(path=path, len_path=len_path)
                    for kid in self._elements
                    for path in kid.all_paths(ko_match)
                }
            )
            return paths

    def abundance(self, ko_match: Dict[str, float]):
        return sum(ko_match.get(ko, 0) for ko in self.list_ko())

    def completeness(self, ko_match: Union[List, Dict, Set]) -> float:
        """Complessness of given match, ko is its dict"""
        count = 0.0
        if self._ko:
            return 1 if self._ko in ko_match else 0
        # multipy elements
        if self._is_chain:
            for element in self._elements:
                count += element.completeness(ko_match)
            return count / len(self._elements)
        # self.is_chain is False
        return max([element.completeness(ko_match) for element in self._elements])


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
