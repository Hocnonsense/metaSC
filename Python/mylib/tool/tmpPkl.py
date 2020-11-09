# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 22:32:22
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-11-09 23:03:12
 * @FilePath: /HScripts/Python/mylib/tool/tmpPkl.py
 * @Description:
    with which will build a tmp pickle file for its function.
"""

import pickle
from functools import wraps


class TmpPkl:
    """ Save a pickle file for the function or block. """

    def __init__(self, PICKLE_FILE_name, force_rewrite=False, **kwargs) -> None:
        self.PICKLE_FILE_name = PICKLE_FILE_name
        self.force_rewrite = force_rewrite
        self.kwargs = kwargs

    def __call__(self, func):
        @wraps(func)
        def wrapTheFunction(*args, **kwargs):
            try:
                if self.force_rewrite:
                    raise FileNotFoundError
                with open(self.PICKLE_FILE_name, "rb") as pi:
                    last_results = pickle.load(pi)
            except (FileNotFoundError, EOFError):
                last_results = func(*args, **kwargs)
                with open(self.PICKLE_FILE_name, "wb") as po:
                    pickle.dump(last_results, po)
            return last_results
        return wrapTheFunction


if __name__ == "__main__":
    @TmpPkl("test.pkl", True)
    def test(i, j):
        return i+j
    print(test(1,3))
