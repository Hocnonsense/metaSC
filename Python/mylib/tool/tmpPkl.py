# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 22:32:22
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-11-14 12:01:24
 * @FilePath: /HScripts/Python/mylib/tool/tmpPkl.py
 * @Description:
    with which will build a tmp pickle file for its function.
        also get some message from pickle: time, path, and desc
        #WARNING# be careful when using @TmpPkl.

 * @TODO:
"""
__version__ = "0.0.3"

import sys
from sys import stderr
import os
import pickle
from functools import wraps
from datetime import datetime
from typing import Callable


class TmpPkl:
    """ Save a pickle file for the function or block. """

    def __init__(self, PICKLE_FILE_name, force_rewrite=False, desc="") -> None:
        self.PICKLE_FILE_name = PICKLE_FILE_name
        self.force_rewrite = force_rewrite
        self.last_results = None
        self.meta = {
            "date": datetime.now(),
            "pwd": sys.argv[0] or os.path.abspath(),
            "desc": desc
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
        return wrappedFunction

    def __enter__(self):
        if not self.force_rewrite:
            try:
                with open(self.PICKLE_FILE_name, "rb") as pi:
                    (self.meta, self.last_results) = pickle.load(pi)
            except (FileNotFoundError, EOFError):
                self.force_rewrite = True
            else:
                print("# load from",
                      os.path.abspath(
                          os.path.expanduser(self.PICKLE_FILE_name)),
                      file=stderr)
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        if self.force_rewrite:
            with open(self.PICKLE_FILE_name, "wb") as po:
                pickle.dump((self.meta, self.last_results), po)
            print("# dump to",
                    os.path.abspath(
                        os.path.expanduser(self.PICKLE_FILE_name)),
                    file=stderr)
        return False

    def __set_desc(self, desc):
        self.meta["desc"] = desc

    desc = property(
        fget=lambda self: self.meta["desc"],
        fset=__set_desc
    )


if __name__ == "__main__":
    @TmpPkl("test.pkl", True)
    def test(i, j):
        """ Test i+j. """
        return i+j
    print(test("yy","sy"))
    with TmpPkl("test.pkl", desc="with i + j") as tmp1:
        tmp1.desc = "with yybsy"
        if tmp1.force_rewrite:
            tmp1.last_results = "ii"+"Vi"
            tmp1.desc["param"] = {"i": "ii", "j": "vi"}
        print(tmp1.meta)
        word = tmp1.last_results
    print(word)
