# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-04 16:36:54
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 14:37:05
 * @FilePath: /HScripts/Python/mylib/tool/func.py
 * @Description: 一个上下文管理器.
 * @TODU:
        1. 为什么 _LOGGER.index() 会使 "\r" 指令无效?
"""


import os
import sys
from datetime import datetime
from contextlib import ContextDecorator
from typing import Any


MAX_LINE_LENGTH = 80


class Func(ContextDecorator):
    """
    New Func for blocking and logging.
    用于块缩进, 段落计时, 输出格式转换, 输出重定向
    @useage
        1.  在需要分块的程序块前添加, 以对该段进行管理
            `with Func():`
        2.  若需要记录时间, 则将该程序块命名为 "blockname", 使用
            `with Func(blockname) as f:
        3.  指定代码块退出条件
            `   f.exit(__name__)`
            等价于
            `if __name__ == "__main__"`
        4.  需要记录输出到与程序块名字相同的文件中, 文件位于 .vscode\\.log 文件夹中
            `   f.open_log()`
        5.  需要获得当前程序块运行时间:
            `   f.time
        6.  需要修改当前最长行长度:
            `   f.set_line_length()`
    """

    def __init__(self, name=None, index=4, **kwargs: Any):
        """
        @param name 标记代码块.  若无则不会缩进
                    输入名字则默认以这个名字作为输出文件的文件名
                    某代码块的名字和用时会额外出现在这个代码块的 log 末尾

        @param index 记录缩进.  默认为 "|   "
                     正整数: 则记录缩进数, 负数为 "|"
                     字符串: 直接作为缩进
                     `0` 或其他: 无缩进
                     `index` 只会作用在窗体输出中, 对文件输出没有作用.

        @param debug 用于输出 bug.  默认为 `True` 即显示错误
                     `False` 或 0: 不显示错误
        """
        self.name = name
        self.log = ""
        self.index = index
        self.no_exception = False
        if "debug" in kwargs:
            self.no_exception = not kwargs["debug"]
        self.start = 0

    def __enter__(self):
        if self.name:
            if self.index and isinstance(self.index, int):
                self.index = "|" + " " * (self.index - 1)
            elif not isinstance(self.index, str):
                self.index = ""
            print("+" + ">" * (len(self.index) - 2) + " {0:15}".format(self.name))
            self.index = _LOGGER.index(self.index)
        self.start = datetime.now()
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        """
        由于格式需要, 结束输出时会回到当前行首.
        可能会丢失部分数据.
        """
        _LOGGER.write_terminal("\r")
        if self.name:
            _LOGGER.index()
            _LOGGER.write_terminal("\r")  # MAGIC! why?
            print("^" + ">" * (len(self.index) - 2) + " {0:15} : {1} seconds".format(self.name, self.time))
        if self.log:
            _LOGGER.close(self.log)
        # 如果为False(None)，异常会被抛出，用户需要进行异常处理。如果为True，则表示忽略该异常。
        return self.no_exception

    @property
    def time(self):
        """Return time spent in this block till now."""
        return (datetime.now() - self.start).total_seconds()

    def exit(self, name, equal="__main__", warning=""):
        """
        Exit block if check is False.
        通过此处退出一定不会报错, 所在行内容清空, 块用时也不会显示
        """
        if not name == equal:
            self.no_exception = True
            if warning:  # 保证自定义警告缩进正确
                print(warning)
            self.name = ""
            _LOGGER.index()
            raise BaseException("٩(๑`н´๑)۶")

    def open_log(self, name=None):
        """
        Open file to write log in.
        每个程序块都可以输出到一个对应的文件中, 但每时刻只能向一个文件中输出.
        返回文件名.
        """
        if self.log:
            _LOGGER.close(self.log)
        self.log = _LOGGER.open(name or self.name or ".log")
        return self.log

    @staticmethod
    def set_line_length(length: int = 80):
        """修改行长度. 若长度为 0, 则同时取消缩进和换行"""
        global MAX_LINE_LENGTH
        MAX_LINE_LENGTH = length


class _Logger(object):
    """LOG 单例模式.  func 模块内部类"""

    def __init__(self, path=".vscode\\.pylog\\"):
        self.path = path
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        self.terminal = sys.stdout
        self.__log = {}
        self.__index = [""]
        sys.stdout = self

    def write(self, message):
        """
        write and format.
        换行会自动添加缩进
        if there is index and set MAX_LINE_LENGTH:
            apply to window stdout
        if log files to write:
            write them without index
        """
        # to window, if index exits
        self.write_terminal(message)
        # to file, if there are
        for log in self.__log.values():
            log.write(message)

    def flush(self):
        """Flush."""
        self.terminal.flush()
        for log in self.__log.values():
            log.flush()

    def index(self, index=None):
        """
        Set indent.
        Return index in/out
        Pop last index if input nothong
        if not index:
            reset sys.stdout
        """
        if index is None:
            return self.__index.pop()
        self.__index.append(index)
        self.terminal.write(index)
        return index

    def write_terminal(self, message):
        """仅向 shell 终端中输出信息"""
        if self.__index and MAX_LINE_LENGTH:
            index = "".join(self.__index)
            lines = str(message).splitlines(True)
            stdout = []
            restlen = MAX_LINE_LENGTH - len(index)
            for line in lines:
                split = len(line) // restlen
                for i in range(split):
                    stdout.append(line[i * restlen:
                                       i * restlen + restlen
                                       ] + "\n" + index)
                lastline = line[split * restlen:] + "\t"
                stdout.append(lastline[:-1])
                if len(lastline.splitlines()) == 2:
                    stdout.append(index)
            self.terminal.write("".join(stdout))
        else:
            self.terminal.write(message)

    def open(self, filename):
        """Open a file to input."""
        if not (filename.startswith("/") or ":" in filename):
            filename = self.path + filename
        self.__log[filename] = self.__log.get(
            filename, open(filename, "a"))
        return filename

    def close(self, filename):
        """Close a file from output."""
        if self.__log:
            log = self.__log.pop(filename)
            log.flush()
            log.close()


_LOGGER = _Logger()


with Func() as f:
    f.exit(__name__)
    with Func("test_1", debug=1) as f1:
        f1.openlog()
        {a: b for a, b in list(map(lambda i: (i, str(i)), range(200000)))}
        print("?\n??")
        with Func("test_2", relocate=False, index=2) as f2:
            f2.openlog("else")
            {a: b for a, b in list(map(lambda i: (i, str(i)), range(200000)))}
            print("???")
            # f2.exit(0, warning="exit 0")
            print("喵喵喵")
        print(
            {a: b for a, b in list(map(lambda i: (i, str(i)), range(20)))},
            end="")
    print("this should not be stored")
