# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-28 21:32:51
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-28 21:46:07
 * @FilePath: /metaSC/PyLib/test/__init__.py
 * @Description:
"""

from pathlib import Path
from typing import Final
import shutil


test_file_path: Final = Path(__file__).parent / "test_temp"
shutil.rmtree(test_file_path, ignore_errors=True)
test_file_path.mkdir(exist_ok=True)
