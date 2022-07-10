# -*- coding: utf-8 -*-
"""
 * @Date: 2020-11-09 22:32:22
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-10 12:18:14
 * @FilePath: /metaSC/PyLib/PyLibTool/tmpPkl.py
 * @Description:
    with which will build a tmp pickle file for its function.
        also get some message from pickle: time, path, and desc
        #WARNING# be careful when using @TmpPkl.
 * @TODO:
"""
__version__ = "0.0.4"

import pickle
from datetime import datetime
from functools import wraps
from pathlib import Path
from sys import argv
from typing import Any, Callable, Dict, Optional, Union

from PyLib.PyLibTool.file_info import md5sum, verbose_import

logger = verbose_import(__name__, __doc__)


class TmpPkl:
    """Save a pickle file for the function or block."""

    __cache: Dict[Callable, Path] = {}  #

    def __init__(
        self,
        PICKLE_FILENAME: Union[str, Path],
        desc: Optional[str] = None,
        force_rewrite=False,
        situ=Path("."),
    ) -> None:
        if situ:
            if not isinstance(situ, Path):
                logger.warning(
                    "var situ is not a Path, may cause problems in future versions"
                )
                situ = Path(situ)
            if situ.is_file():
                pickle_filename = situ.parent / PICKLE_FILENAME
            else:
                pickle_filename = situ / PICKLE_FILENAME
        else:
            pickle_filename = Path(PICKLE_FILENAME)
        self.PICKLE_FILENAME = pickle_filename.expanduser().absolute()
        self.force_rewrite = force_rewrite
        self.__force_rewrite = force_rewrite
        self.last_results: Any = None
        self.meta = {
            "date": datetime.now(),
            "pwd": Path(argv[0] or ".").absolute(),
            "desc": desc or "",
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
                logger.info(f"# load from {self.PICKLE_FILENAME} ... ")
                with open(self.PICKLE_FILENAME, "rb") as pi:
                    (self.meta, self.last_results) = pickle.load(pi)
            except (FileNotFoundError, EOFError):
                self.force_rewrite = True
                logger.critical(f"\r# load from {self.PICKLE_FILENAME} ... failed")
            else:
                logger.critical(
                    f"\r# load from {self.PICKLE_FILENAME} ... finished, "
                    f"md5sum {md5sum(file=self.PICKLE_FILENAME)}"
                )
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        if self.force_rewrite and self.last_results:
            logger.info(f"# dump to {self.PICKLE_FILENAME} ... ")
            with open(self.PICKLE_FILENAME, "wb") as po:
                pickle.dump((self.meta, self.last_results), po)
            logger.critical(
                f"\r# dump to {self.PICKLE_FILENAME} ... finished, "
                f"md5sum {md5sum(file=self.PICKLE_FILENAME)}"
            )
        self.force_rewrite = self.__force_rewrite
        return False

    def __set_desc(self, desc):
        self.meta["desc"] = desc

    desc = property(fget=lambda self: self.meta["desc"], fset=__set_desc)

    # unused
    # @classmethod
    # def __get_cache(cls):
    #    return cls.__cache

    def __show_cache(self):
        return {k: v.is_file() for k, v in self.__class__.__cache.items()}

    cache = property(fget=__show_cache)
