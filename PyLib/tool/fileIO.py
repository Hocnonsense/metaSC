# -*- coding: utf-8 -*-
"""
 * @Date: 2020-08-20 00:16:45
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 14:27:10
 * @FilePath: /HScripts/Python/mylib/tool/fileIO.py
 * @Description: 文件输入输出
"""

import os
import io
import json
from typing import Any, AnyStr


def write_file(file_name: str, content: AnyStr) -> None:
    """Write string `content` as a text file."""
    with open(file_name, "w") as f:
        f.write(content)


def write_json_file(file_name: str, obj: Any, indent: int = 4) -> None:
    """Write `obj` to a file in a pretty JSON format. This supports Unicode."""
    # Indentation puts a newline after each ',' so suppress the space there. | "," 后跟 "\n", 所以 ", " 没有必要
    message = json.dumps(obj, ensure_ascii=False, indent=indent,
                         separators=(',', ': '), sort_keys=True) + '\n'
    write_file(file_name, message.encode('utf-8'))


def read_json_file(file_name: str) -> Any:
    """Read and parse JSON file. This supports Unicode."""
    with io.open(file_name, encoding='utf-8') as f:
        return json.load(f)


def is_exist(file_name: str) -> None:
    """
    判断文件是否存在, 提醒是否需要覆盖写入文件.
    """
    if os.path.exists(file_name):
        print(file_name + " 已存在, 是否覆盖? [sty]")
        mode = input()
        if mode.lower()[0] in "sty":
            with open(file_name, 'w') as clear:
                clear.write("")
        elif mode.lower()[0] in "bfn":
            print("是否将已存在的文件重命名? [sty]/[bfn]")
            mode = input()
            if mode.lower()[0] in "sty":
                i = 1
                while os.path.exists(file_name + "(" + str(i) + ")"):
                    i += 1
                os.rename(file_name, file_name + "(" + str(i) + ")")
            elif mode.lower()[0] in "bfn":
                print("Warning: 已存在同名文件, 内容未知. ")
            else:
                raise FileExistsError("那你想要怎样嘛")
        else:
            raise FileExistsError("写文件进程退出")
    else:
        print("成功创建" + file_name)


def eopen(file_name: str, *args, **kwargs):
    """ @descrption如果是重写, 则提醒是否需要覆盖原文件.  "w", "wb" 会被提醒, 但 "wt" 不会
    """
    if args and args[0] in ["w", "wb", ]:
        is_exist(file_name)
    return open(file_name, *args, **kwargs)


eopen.__doc__ += open.__doc__
