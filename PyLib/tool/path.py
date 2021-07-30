# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-05 11:03:35
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-07-27 15:54:17
 * @FilePath: /metaSC/PyLib/tool/path.py
 * @Description: 文件和文件名操作子
 * @TODO:
"""


import errno
import os
import re
from typing import Iterable, List, Tuple


"""
filepath.py
File and file_name path utilities.
"""


def makedirs(path: str, *paths: Iterable[str]) -> str:
    """ @description 创建路径
            Join one or more path components, make that directory path (using the
            default mode 0o0777), and return the full path.
        @return 完整路径
            Raise OSError if it can't achieve the result (e.g. the containing directory
            is readonly or the path contains a file); not if the directory already
            exists.
    """
    full_path = os.path.join(path, *paths)

    os.makedirs(full_path, exist_ok=True)

    return full_path


def verify_file_exists(file_path: str, message: str = '') -> None:
    """ @description assert 文件存在
            Raise an FileNotFoundError if file_path isn't an existing file."""
    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            errno.ENOENT, 'Missing file "{}".  {}'.format(file_path, message))


def verify_dir_exists(dir_path: str, message: str = '') -> None:
    """ @description assert 路径存在
            Raise an NotADirectoryError if dir_path isn't an existing directory."""
    if not os.path.isdir(dir_path):
        raise NotADirectoryError(
            errno.ENOENT, 'Missing dir "{}".  {}'.format(dir_path, message))


def list_dir(dir: str, filter: Tuple) -> List:
    """ List file or direcory in `dir`. """
    dir_list = []
    if isinstance(filter, str):
        filter = re.compile(filter)
    for file in sorted(os.listdir(dir)):
        if re.match(filter, file):
            dir_list.append(file)
    return dir_list


def fine_all_files(dir: str, filter) -> List:
    """ List all files. """
    # : re.Pattern is unavalible in python3.6
    files_list = []
    for parent, dirs, files in os.walk(dir):
        file_list = []
        for file in files:
            if re.match(filter, file):
                file_list.append(file)
        files_list.append((parent, file_list))
    return files_list


def abs_user_path(dir: str) -> str:
    """ Recognize `~/` and modify it. """
    return os.path.abspath(os.path.expanduser(dir))


class Path:
    """
    这个类是帮助管理路径的.
    最初用在 pdfhelper 这个项目中.
    输入是文件名, 路径, 存储文件夹的名字, 默认路径为当前项目, 存储文件名与文件名相同

    __call__ 直接调用 eopen, 就是 open 的升级版, 能提示文件存在情况, 方便决定是否替换或者别的.

    now_path: 当前运行的主 Python 文件所在的目录
    source_name: 当前路径管理的名字
    source_path: 当前的路径
    """
    __dict = {}
    __rootpath = os.path.dirname(__file__)

    def __init__(self, source_path, desc):
        """
        @param source_path 当前路径, 若无, 则为调用程序根目录.
        @param desc 当前路径管理的描述.  推荐为要读取的文件的名字
        """
        if source_path in Path.__dict:
            raise AttributeError("source path already exits!")
        if ":" not in source_path and source_path[0] != "/":
            raise NotADirectoryError("only accept absolute directory!")
        self.__sourcepath = source_path
        self.__desc = desc

    @classmethod
    def cd(cls, *source_path, desc=""):
        # param source_path check
        if not source_path:
            source_path = cls.root_path
        source_path = os.path.join(*source_path)

        if os.path.isfile(source_path):
            source_path = os.path.dirname(source_path)
        elif not os.path.isdir(source_path):
            source_path = makedirs(source_path)

        source_path = os.path.abspath(source_path)

        # source_path defined check
        if source_path in cls.__dict:
            if desc:
                print("Warning: an Path already exits!")
            return cls.__dict[source_path]

        path = cls(source_path, desc)
        cls.__dict[source_path] = path

        return path

    @classmethod
    @property
    def root_path(cls):
        """调用文件的根路径"""
        return cls.__rootpath

    @property
    def source_path(self):
        """
        当前路径, 若无, 则为调用程序所在位置.
        只读
        """
        return self.__sourcepath

    @property
    def desc(self):
        """
        当前路径管理的名字.
        只读
        """
        return self.__desc


"""
    规定:
        路径用 "\\" 分隔, 以 "\\"结束
"""
if __name__ == "__main__":
    # pylint: disable = invalid-name
    source_path = "C:\\Users\\Hwrn"
    desc = ".condarc"
    file1 = Path.cd(source_path, desc=desc)
    print(file1.source_path, file1.desc)
    print(file1.root_path, Path.root_path)
    try:
        Path.cd(source_path, desc="another")
        Path(source_path, desc="fault")
    except AttributeError as e:
        print(e)
    try:
        Path("mylib", "fault1")
    except NotADirectoryError as e:
        print(e)
    file2 = Path.cd("mylib", desc="__init__.py")
    print(file2.source_path, file2.desc)
