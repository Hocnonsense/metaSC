from _typeshed import Incomplete as Incomplete
from typing import Optional

logger: Incomplete

class Timer:
    timer_start: Incomplete
    initial_checkpoint_key: Incomplete
    last_checkpoint_key: Incomplete
    checkpoints: Incomplete
    num_checkpoints: int
    required_completion_score: Incomplete
    score: Incomplete
    complete: bool
    last_eta: Incomplete
    last_eta_timestamp: Incomplete
    scores: Incomplete
    def __init__(self, required_completion_score: Incomplete | None = ..., initial_checkpoint_key: int = ..., score: int = ...) -> None: ...
    def timestamp(self) -> None: ...
    def timedelta_to_checkpoint(self, timestamp, checkpoint_key: Incomplete | None = ...): ...
    def make_checkpoint(self, checkpoint_key: Incomplete | None = ..., increment_to: Incomplete | None = ...): ...
    def gen_report(self, title: str = ...) -> None: ...
    def gen_dataframe_report(self) -> None: ...
    def gen_file_report(self, filepath) -> None: ...
    def calculate_time_remaining(self, infinite_default: str = ...): ...
    def eta(self, fmt: Incomplete | None = ..., zero_padding: int = ...): ...
    def time_elapsed(self, fmt: Incomplete | None = ...): ...
    def format_time(self, timedelta, fmt: Optional[str] = ..., zero_padding: int = ...): ...

class TimeCode:
    single_line_prefixes: Incomplete
    quiet: Incomplete
    suppress_first: Incomplete
    s_msg: Incomplete
    f_msg: Incomplete
    def __init__(self, success_msg: Incomplete | None = ..., sc: str = ..., fc: str = ..., failure_msg: Incomplete | None = ..., quiet: bool = ..., suppress_first: int = ...) -> None: ...
    timer: Incomplete
    def __enter__(self) -> None: ...
    time: Incomplete
    def __exit__(self, exception_type, exception_value, traceback) -> None: ...
