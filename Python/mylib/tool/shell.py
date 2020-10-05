# -*- coding: utf-8 -*-
"""
 * @Date: 2020-08-19 18:13:55
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 15:29:51
 * @FilePath: /HScripts/Python/mylib/tool/shell.py
 * @Description:
"""

import os
import subprocess
import sys
from typing import Optional, Sequence, Tuple


def printsh(value, flow=0, end="\n"):
    """ Help to write in file """
    value = value + end
    if flow != 2:  # write to stdout
        sys.stdout.write(value)
    if flow != 1:  # write to stderr
        sys.stderr.write(value)


def runsh(tokens: Sequence[str], trim: bool = True) -> str:
    """ @description: 在管道中运行 shell 命令
            Run a shell command-line (in token list form) and return its output.
            This does not expand filename patterns or environment variables or do other
            shell processing steps.

            This sets environment variables `PATH` and (if available) `LD_LIBRARY_PATH`.
            Sherlock needs the latter to find libcrypto.so to run `git`.

        @param
            tokens: The command line as a list of string tokens.
            trim: Whether to trim off trailing whitespace. This is useful
                because the subprocess output usually ends with a newline.
        @returns 命令行输出
            The command's output string.
    """
    if ">" in "".join(tokens):
        raise NotImplementedError("runsh is using subprecess.Popen, cannot recognize popens!")

    environ = {
        "PATH": os.environ["PATH"],
        # 临时设置动态链接库的地址
        "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
        "PYTHONHASHSEED": "1"
    }
    out = subprocess.Popen(tokens, stdout=subprocess.PIPE,
                           env=environ)
    out = out.communicate()[0]
    if trim:
        out = out.rstrip()
    return out


def runsh_safe(tokens: Tuple[str, list], trim: bool = True) -> Optional[str]:
    """Run a shell command-line string and return its output. This does not
    expand filename patterns or environment variables or do other shell
    processing steps.

    This sets environment variables `PATH` and (if available) `LD_LIBRARY_PATH`.
    Sherlock needs the latter to find libcrypto.so to run `git`.

    Args:
        line: The command line as a string.
        trim: Whether to trim off trailing whitespace. This is useful
            because the subprocess output usually ends with a newline.
    Returns:
        The command's output string, or None if it couldn't even run.
    """
    if isinstance(tokens, str):
        tokens = tokens.split()
    if sys.platform == 'win32':
        tokens = ["cmd", "/C"] + tokens
    try:
        return runsh(tokens, trim=trim)
    # pylint: disable = broad-except
    except Exception as e:
        print('failed to run command line {}: {}'.format(tokens, e))
        return None
