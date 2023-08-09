# -*- coding: utf-8 -*-
"""
 * @Date: 2022-10-04 19:48:03
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-09 15:30:03
 * @FilePath: /metaSC/PyLib/PyLibTool/source.py
 * @Description:
    source file by given name.
"""

import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType


def source(file: Path) -> ModuleType:
    """
    Loads, executes, then returns the module at input file.
    Raises:
        IOError: If the module cannot be loaded. This
            exception is chosen.
    Returns:
        types.ModuleType: The initialised module.
    Reference:
        https://github.com/wntrblm/nox/pull/498/files
    """
    # make sure file is a Path object
    file = Path(file).expanduser().absolute()

    module_name = file.name.rsplit(".", 1)[0]
    spec = importlib.util.spec_from_file_location(module_name, file)
    if not spec:
        raise IOError(f"Could not get module spec from {file}")

    module = importlib.util.module_from_spec(spec)
    if not module:
        raise IOError(f"{module_name} is not a valid python module.")

    sys.modules[module_name] = module
    loader = spec.loader
    if not loader:  # pragma: no cover
        raise IOError(f"Could not get module loader for {module_name}")

    loader.exec_module(module)  # type: ignore
    return module


def source_env(env_var: str):
    env_value = os.environ.get(env_var, "")
    if not env_value:
        return
    file = Path(env_value).expanduser().absolute()
    if not file.exists():
        raise FileNotFoundError(f"please check file {file}")
    return source(file)
