# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-30 20:53:37
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 18:50:45
 * @FilePath: /metaSC/PyLib/tool/wildcards/common.py
 * @Description:
"""

import collections.abc
import inspect
from pathlib import Path
from types import FrameType
from typing import Any, Optional, Union, overload

PathLike = Union[str, Path]


class Unformattable(ValueError):
    def __init__(self, errormsg="This cannot be used for formatting"):
        self.errormsg = errormsg


class Pattern:
    @overload
    def __init__(self, pattern: "Pattern") -> None:
        pass

    @overload
    def __init__(self, pattern: PathLike) -> None:
        pass

    @overload
    def __init__(self, pattern: str, _type: type) -> None:
        pass

    def __init__(self, pattern, _type: Optional[type] = None) -> None:
        if isinstance(pattern, Pattern):
            self.pattern: str = pattern.pattern
            self._type: type = pattern._type
        else:
            self.pattern = str(pattern)
            if _type is not None:
                self._type = _type
            elif isinstance(pattern, Path) or isinstance(pattern, str):
                self._type = type(pattern)
            else:
                self._type = str

    def __repr__(self) -> str:
        return f"{type(self).__name__}[{self._type}]({str(self)})"

    def __str__(self) -> str:
        return self.pattern

    @property
    def real(self):
        return self._type(self.pattern)

    @property
    def names(self):
        names: list[str] = []
        return names

    def __eq__(self, __value: object) -> bool:
        if isinstance(__value, Pattern):
            return self.pattern == __value.pattern and isinstance(
                __value._type, self._type
            )
        else:
            return self.pattern == str(__value) and isinstance(__value, self._type)


def stepout(nstep=1):
    if nstep == 0:
        return {}
    frame: FrameType = inspect.currentframe().f_back  # type: ignore
    while nstep > 1:
        print("_stepout", frame.f_code.co_filename, flush=True)
        if not frame.f_back:
            break
        frame = frame.f_back
        nstep -= 1

    # add local variables from calling rule/function
    return {**frame.f_globals, **frame.f_locals}


def flatten(wildcards: dict):
    for wildcard, values in wildcards.items():
        if isinstance(values, str) or not isinstance(values, collections.abc.Iterable):
            values = [values]
        yield [(wildcard, value) for value in values]


def extract_args_kwargs(comb: dict[Union[str, int], Any]):
    args: list[Optional[Any]] = []
    kwargs: dict[str, Any] = {}
    for k, v in comb.items():
        if isinstance(k, int):
            if k >= len(args):
                args += (k - len(args) + 1) * [None]
            args[k] = v
            continue
        kwargs[k] = v
    return args, kwargs
