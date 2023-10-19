from typing import Callable, Dict, List, TextIO, Union

def load_KEGG_module_raw(source: Union[str, TextIO], cache_path: str = ...) -> Dict[str, List[str]]: ...
def read_brite_json(brite_doc) -> None: ...
def load_brite(source: Union[str, TextIO], cache_path: str = ...): ...
def cached(func: Callable[[str], TextIO], source: str, filename: Union[None, str] = ...): ...
