from _typeshed import Incomplete
from enum import Enum

__maintainer__: str
__status__: str

class ShellColors(Enum):
    gray: str
    red: str
    green: str
    yellow: str
    blue: str
    magenta: str
    cyan: str
    white: str
    crimson: str

class ShellWeights(Enum):
    normal: str
    bold: str

tty_colors: Incomplete

def color_text(text, color: str = ..., weight: str = ...): ...
