from PyLib.PyLibTool.file_info import verbose_import as verbose_import
from _typeshed import Incomplete
from pathlib import Path
from typing import TextIO

logger: Incomplete

def prodigal(genome: Path, out_dirname: Path = ..., out_basename: str = ..., filter: bool = ..., meta: bool = ..., min_faa_len: int = ...): ...
def prodigal_filter(out_prefix: Path, filter_prefix: Path = ..., min_faa_len: int = ...): ...
def prodigal_faa2gff_iter(pfain: TextIO): ...
def prodigal_faa2gff_cmd(gene_prefix: Path): ...
def run(loglevel: str, genome: Path, gene_prefix: str, filter: bool = ..., meta: bool = ..., min_faa_len: int = ...): ...
