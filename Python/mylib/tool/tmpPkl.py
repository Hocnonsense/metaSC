# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 22:32:22
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-11-11 21:05:58
 * @FilePath: /HScripts/Python/mylib/tool/tmpPkl.py
 * @Description:
    with which will build a tmp pickle file for its function.
 * @TODO: also get some message from pickle
"""
__version__ = "0.0.2"

from sys import stderr
import os
import pickle
from functools import wraps


class TmpPkl:
    """ Save a pickle file for the function or block. """

    def __init__(self, PICKLE_FILE_name, force_rewrite=False, **kwargs) -> None:
        self.PICKLE_FILE_name = PICKLE_FILE_name
        self.force_rewrite = force_rewrite
        self.kwargs = kwargs
        self.last_results = None

    def __call__(self, func):
        @wraps(func)
        def wrappedFunction(*args, **kwargs):
            with self as tmppkl:
                if tmppkl.force_rewrite:
                    tmppkl.last_results = func(*args, **kwargs)
            return tmppkl.last_results
        return wrappedFunction

    def __enter__(self):
        if not self.force_rewrite:
            try:
                with open(self.PICKLE_FILE_name, "rb") as pi:
                    self.last_results = pickle.load(pi)
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
                pickle.dump(self.last_results, po)
            print("# dump to",
                    os.path.abspath(
                        os.path.expanduser(self.PICKLE_FILE_name)),
                    file=stderr)
        return False


if __name__ == "__main__":
    @TmpPkl("test.pkl", )
    def test(i, j):
        return i+j
    print(test("yy","sy"))
    with TmpPkl("test.pkl") as tmp1:
        if tmp1.force_rewrite:
            tmp1.last_results = "yy"+"bsy"
        word = tmp1.last_results
    print(word)
