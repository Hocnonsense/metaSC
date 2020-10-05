# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-06 16:18:59
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 15:31:33
 * @FilePath: /HScripts/Python/mylib/tool/Logger.py
 * @Description:
"""

import sys
import time


class Logger(object):
    """
    *args, **kwargs for open(filename, ... )
    sys.stdout = Logger( _yourlogfilename_ )

    then everything in python will save in your log file.
"""
    def __init__(self, *args, filename="Default.log", **kwargs):
        self.terminal = sys.stdout
        self.log = open(filename, *args, **kwargs)

    def write(self, message):
        # message= message.encode('unicode_escape').decode('utf-8')
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass


defaultTitle = time.strftime(
    "new LOG at %Y-%m-%d %H:%M:%S :", time.localtime())


def savelog(filename="programLog.txt", title=defaultTitle):
    sys.stdout = Logger(filename, "a")
    print(title)
