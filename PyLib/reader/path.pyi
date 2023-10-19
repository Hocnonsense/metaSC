from _typeshed import Incomplete

class path:
    default_file: str
    path: Incomplete
    def __init__(self, sourse_path: str = ..., **kwargs) -> None: ...
    def get(self, filename: str = ..., *filenames) -> str: ...
