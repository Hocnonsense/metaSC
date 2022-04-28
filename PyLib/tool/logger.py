# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-06 16:18:59
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-28 23:01:47
 * @FilePath: /metaSC/PyLib/tool/logger.py
 * @Description:
"""

import os
import sys


class Tee:
    """
    everything in python will save in your log file.
    *args, **kwargs for open(filename, ... )
    >>> sys.stdout = Tee(_yourlogfilename_)
    """

    def __init__(self, *args, filename="Default.log", **kwargs):
        self.log = open(filename, *args, **kwargs)
        self.logger=sys.stdout

    def write(self, message):
        # message= message.encode('unicode_escape').decode('utf-8')
        self.logger.write(message)
        self.log.write(message)

    def flush(self):
        self.logger.flush()
        self.log.flush()


class Redirector:
    def __init__(self, stderr=os.devnull, stdout=os.devnull):
        self.stdout = stdout
        self.stderr = stderr

    def __enter__(self):
        sys.stderr.flush()
        self.old_stderr = sys.stderr
        if isinstance(self.stderr, str):
            sys.stderr = open(self.stderr, "w")
        else:
            sys.stderr = self.stderr

        sys.stdout.flush()
        self.old_stdout = sys.stdout
        if isinstance(self.stdout, str):
            sys.stdout = open(self.stdout, "w")
        else:
            sys.stdout = self.stdout

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.flush()
        sys.stderr = self.old_stderr
        sys.stdout.flush()
        sys.stdout = self.old_stdout
