# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-03 13:09:46
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-03 15:17:04
 * @FilePath: /MyScripts/Python/mylib/tool/parseArgs.py
 * @Description:
        Common code for scripts that manually run simulation and analysis operations
        from the command line (e.g. outside of Fireworks workflows). This has code for
        finding the simulation dir, variant subdirs, etc.

        Run with '-h' for command line help.
        Set PYTHONPATH when running this.
"""

# 详见 https://docs.python.org/zh-cn/3/library/argparse.html
import argparse
from typing import Any


def str_to_bool(s: str) -> bool:
    """
        @description: 识别 "是否" 类型的参数
        @param s:
            true, false, 1, or 0.
        @return {type}
    """
    s = s.lower()
    if s not in {'true', 'false', '1', '0'}:
        raise ValueError('Expected a bool, not %s' % s)
    return s in {'true', '1'}


def define_parameters(parser: argparse.ArgumentParser) -> None:
    """@description: 参数设置
        Define command line parameters. This base method defines a --verbose
        flag. Overrides should call super.

        Examples include positional arguments
            `parser.add_argument('variant', nargs='?',
            help='Simulation variant.')`
        options
            `parser.add_argument('--seed', default='000000',
            help='Simulation seed.')`.
        and flags
            `parser.add_argument('--verbose', action='store_true',
            help='Enable verbose logging.')`.
    """
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose logging.')


def define_option(parser: argparse.ArgumentParser,
                  name: str,
                  default: Any,
                  help: str,
                  short_name: str = ""
                  ) -> None:
    """@description: 可选参数设置
        Add an option with the given name and datatype to the parser.
    """
    datatype = type(default)
    if short_name:
        parser.add_argument('-' + short_name, '--' + name,
                            type=type(default), default=default,
                            help='({}; {}) {}'.format(
                                datatype.__name__, default, help)
                            )
    else:
        parser.add_argument('--' + name,
                            type=type(default), default=default,
                            help='({}; {}) {}'.format(
                                datatype.__name__, default, help)
                            )


def define_parameter_bool(parser: argparse.ArgumentParser,
                          name: str,
                          default: Any,
                          help: str
                          ) -> None:
    """@description: "bool" 类型的参数设置. 有默认值
        Add a boolean option parameter to the parser. The CLI input can be
        `--name`, `--no_name`, `--name true`, `--name false`, `--name 1`,
        `--name 0`, `--name=true`, etc. The default can be True or False, and
        changing it won't affect any of those explicit input forms. This method
        adds the default value to the help text.
    """
    default = bool(default)
    examples = 'true or 1' if default else 'false or 0'
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--' + name,
                       nargs='?', default=default, const='true',  # needed for nargs='?'
                       type=str_to_bool,
                       help='({}; {}) {}'.format('bool', examples, help))
    group.add_argument('--no_' + name, dest=name, action='store_false',
                       help='Like {}=0'.format(name))


def send_parser(description):
    """ Send an argparse.ArgumentParser to user.
       @param description: description
        """
    return argparse.ArgumentParser(description=description)
