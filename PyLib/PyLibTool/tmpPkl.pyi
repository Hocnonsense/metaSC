from PyLib.PyLibTool.file_info import md5sum as md5sum, verbose_import as verbose_import
from _typeshed import Incomplete
from pathlib import Path
from typing import Callable, Optional, Union

logger: Incomplete

class TmpPkl:
    PICKLE_FILENAME: Incomplete
    force_rewrite: Incomplete
    last_results: Incomplete
    meta: Incomplete
    def __init__(self, PICKLE_FILENAME: Union[str, Path], desc: Optional[str] = ..., force_rewrite: bool = ..., situ=...) -> None: ...
    def __call__(self, func: Callable): ...
    def __enter__(self): ...
    def __exit__(self, exc_type: Incomplete | None = ..., exc_val: Incomplete | None = ..., exc_tb: Incomplete | None = ...): ...
    desc: Incomplete
    @classmethod
    def clean(cls, tmppkl_or_function: Union['TmpPkl', Callable]): ...
    cache: Incomplete
