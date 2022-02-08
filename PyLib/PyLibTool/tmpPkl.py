# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 22:32:22
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-08 17:36:38
 * @FilePath: /metaSC/PyLib/PyLibTool/tmpPkl.py
 * @Description:
    with which will build a tmp pickle file for its function.
        also get some message from pickle: time, path, and desc
        #WARNING# be careful when using @TmpPkl.
 * @TODO:
"""
__version__ = "0.0.4"

import sys
from sys import stderr
import os
import pickle
from functools import wraps
from datetime import datetime
from typing import Callable, Dict, Tuple


class TmpPkl:
    """Save a pickle file for the function or block."""

    __cache: Dict[Callable, str] = {}  #

    def __init__(self, PICKLE_FILENAME, desc="", force_rewrite=False, situ="") -> None:
        if situ:
            if os.path.isfile(situ):
                situ = os.path.dirname(situ)
            PICKLE_FILENAME = os.path.join(situ, PICKLE_FILENAME)
        else:
            PICKLE_FILENAME = os.path.expanduser(PICKLE_FILENAME)
        self.PICKLE_FILENAME = os.path.abspath(PICKLE_FILENAME)
        self.force_rewrite = force_rewrite
        self.last_results = None
        self.meta = {
            "date": datetime.now(),
            "pwd": sys.argv[0] or os.getcwd(),
            "desc": desc,
        }

    def __call__(self, func: Callable):
        """
        * @description: It is studpid, be careful when using this feature.
        * @param {*} self
        * @param {Callable} func
        * @return {*}
        """

        @wraps(func)
        def wrappedFunction(*args, **kwargs):
            with self as tmppkl:
                if tmppkl.force_rewrite:
                    tmppkl.last_results = func(*args, **kwargs)
            return tmppkl.last_results

        self.meta["desc"] = func.__doc__

        self.__class__.__cache[wrappedFunction] = self.PICKLE_FILENAME

        return wrappedFunction

    def __enter__(self):
        if not self.force_rewrite:
            try:
                print(f"# load from {self.PICKLE_FILENAME} ... ", end="", file=stderr)
                with open(self.PICKLE_FILENAME, "rb") as pi:
                    (self.meta, self.last_results) = pickle.load(pi)
            except (FileNotFoundError, EOFError):
                self.force_rewrite = True
                print(f"failed", file=stderr)
            else:
                print(f"finished", file=stderr)
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        if self.force_rewrite and self.last_results:
            print(f"# dump to {self.PICKLE_FILENAME} ... ", end="", file=stderr)
            with open(self.PICKLE_FILENAME, "wb") as po:
                pickle.dump((self.meta, self.last_results), po)
            print(f"finished", file=stderr)
        return False

    def __set_desc(self, desc):
        self.meta["desc"] = desc

    desc = property(fget=lambda self: self.meta["desc"], fset=__set_desc)

    @classmethod
    def __get_cache(cls):
        return cls.__cache

    def __show_cache(self):
        return {k: os.path.isfile(v) for k, v in self.__class__.__cache.items()}

    cache = property(fget=__show_cache)


if __name__ == "__main__":

    @TmpPkl("test.pkl", True)
    def test(i, j):
        """Test i+j."""
        return i + j

    print(test("yy", "sy"))
    with TmpPkl("test.pkl", desc="with i + j") as tmp1:
        tmp1.desc = "with yybsy"
        if tmp1.force_rewrite:
            tmp1.last_results = "ii" + "Vi"
            tmp1.desc = {"i": "ii", "j": "vi"}
        print(tmp1.meta)
        word = tmp1.last_results
    print(tmp1.cache[test], word)
