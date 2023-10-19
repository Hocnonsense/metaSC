from _typeshed import Incomplete as Incomplete
from pathlib import Path
from typing import Any, Union, overload

PathLike = Union[str, Path]

class Unformattable(ValueError):
    errormsg: Incomplete
    def __init__(self, errormsg: str = ...) -> None: ...

class Pattern:
    @overload
    def __init__(self, pattern: Pattern) -> None: ...
    @overload
    def __init__(self, pattern: PathLike) -> None: ...
    @overload
    def __init__(self, pattern: str, _type: type) -> None: ...
    @property
    def real(self) -> None: ...
    @property
    def names(self) -> None: ...
    def __eq__(self, __value: object) -> bool: ...

def stepout(nstep: int = ...): ...
def flatten(wildcards: dict): ...
def extract_args_kwargs(comb: dict[Union[str, int], Any]): ...
