# -*- coding: utf-8 -*-
"""
 * @Date: 2021-04-26 13:23:04
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-04-26 15:17:24
 * @FilePath: /metaSC/PyLib/reader/path.py
 * @Description:
"""


import os


class path:
    def __init__(self, sourse_path='', **kwargs) -> None:
        self.default_file = ''
        if not sourse_path:
            sourse_path = os.path.abspath('')
        else:
            if sourse_path.startswith('~'):
                sourse_path = os.path.expanduser(sourse_path)
            if os.path.isdir(sourse_path):
                self.path = sourse_path
            elif os.path.isfile(sourse_path):
                self.default_file = os.path.basename(sourse_path)
                sourse_path = os.path.dirname(sourse_path)
            else:
                raise FileNotFoundError
        self.path = sourse_path

    def get(self, filename='', *filenames) -> str:
        """
        * @description:
        * @param {*} filename
        * @return {*}
        """
        if filename == '':
            filename = self.default_file
        if not filename:
            raise FileNotFoundError
        return os.path.abspath(os.path.join(self.path, filename, *filenames))
