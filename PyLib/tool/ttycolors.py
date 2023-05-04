# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-28 23:23:30
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-04 10:55:21
 * @FilePath: /metaSC/PyLib/tool/ttycolors.py
 * @Description:
    TTY colors
"""


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"

import sys
from enum import Enum


class ShellColors(Enum):
    gray = "30"
    red = "31"
    green = "32"
    yellow = "33"
    blue = "34"
    magenta = "35"
    cyan = "36"
    white = "37"
    crimson = "38"


class ShellWeights(Enum):
    normal = "0"
    bold = "1"


tty_colors = {
    i.name: {j.name: f"\033[{j.value};{i.value}m%s\033[0m" for j in ShellWeights}
    for i in ShellColors
}


def color_text(text, color="crimson", weight="bold"):
    if sys.stdout.isatty():
        return tty_colors[color][weight] % text
    else:
        return text
