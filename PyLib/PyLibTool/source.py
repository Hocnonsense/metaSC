# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-04 19:48:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-04 19:56:06
 * @FilePath: /metaSC/PyLib/PyLibTool/source.py
 * @Description:
    source file by given name.
"""

import os
from pathlib import Path
from importlib.machinery import SourceFileLoader


def source(file: Path):
    # make sure file is a Path object
    file = Path(file).expanduser().absolute()

    module_name = file.name.rsplit(".", 1)[0]
    module = SourceFileLoader(module_name, str(file)).load_module()
    return module


def source_env(env_var: str):
    env_value = os.environ.get(env_var, "")
    if not env_value:
        return
    file = Path(env_value).expanduser().absolute()
    if not file.exists():
        raise FileNotFoundError(f"please check file {file}")
    return source(file)
