# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-03 13:09:46
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-03 15:32:28
 * @FilePath: /HScripts/Python/mylib/tool/scriptBase.py
 * @Description:
        Common code for scripts that manually run simulation and analysis operations
        from the command line (e.g. outside of Fireworks workflows). This has code for
        finding the simulation dir, variant subdirs, etc.

        Run with '-h' for command line help.
        Set PYTHONPATH when running this.
"""

# 继承了 abc.ABCMeta 的类, 可以用 @abc.abstractmethod 方法创建虚拟方法. https://www.jianshu.com/p/06c960020322
import abc
# 详见 https://docs.python.org/zh-cn/3/library/argparse.html
import argparse
import os
import datetime
import pprint as pp
import time

try:
    from . import parseArgs as pA
except ImportError:
    # pyright: reportMissingImports=false
    import parseArgs as pA


class ScriptBase(metaclass=abc.ABCMeta):
    """
        @description: Abstract base class for scripts. | 脚本基类
            `description()` describes the script,
            `define_parameters()` defines its command line parameters,
                `define_option` and `define_parameter_bool` can be used to define parameters
                `"verbose"=true` is the only default param
            `parse_args()` parses the command line args,
            `run()` does the work,
            `cli()` is the driving Command-Line Interpreter.
    """
    def description(self):
        """Describe the command line program. This defaults to the class name."""
        return type(self).__name__

    # parse args
    def define_parameters(self, *args, **kwargs):
        return pA.define_parameters(self.parser, *args, **kwargs)

    def modify_args(self, args):
        return args

    @abc.abstractmethod
    def run(self, args):
        # type: (argparse.Namespace) -> None
        """Run the operation with the given arguments. If args.verbose,
        overrides can do verbose logging.
        """
        raise NotImplementedError("ScriptBase subclass must implement run()")

    # not recommand to change
    def define_option(self, *args, **kwargs):
        """ name: str, default: Any, help: str, short_name: str = "" """
        return pA.define_option(self.parser, *args, **kwargs)

    def define_parameter_bool(self, *args, **kwargs):
        """ name: str, default: Any, help: str """
        return pA.define_parameter_bool(self.parser, *args, **kwargs)

    def help(self):
        """Return help text for the Command Line Interface. This defaults to a
            string constructed around `self.description()`.
        """
        return 'Run {}.'.format(self.description())

    def parse_args(self) -> argparse.Namespace:
        """Parse the command line args: Construct an ArgumentParser, call
        `define_parameters()` to define parameters including subclass-specific
        parameters, use it to parse the command line into an
        `argparse.Namespace`, and return that.

        When overriding, first call super().

        (A `Namespace` is an object with attributes and some methods like
        `__repr__()` and `__eq__()`. Call `vars(args)` to turn it into a dict.)
        """
        self.parser = pA.send_parser(self.help())
        self.define_parameters()
        args = self.parser.parse_args()

        return args

    def cli(self):
        """Command Line Interpreter: parse_args() then run(). This also prints
        a starting message (including args.sim_path if defined) and an ending
        message (including the elapsed run time).
        """
        start_wall_sec = time.time()
        print('{} > {} > {}'.format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_wall_sec)),
            os.getcwd(),
            self.description()))

        args = self.parse_args()
        pp.pprint({'Arguments': vars(args)})

        start_process_sec = time.perf_counter()
        self.run(args)  # !RUNNING! #
        elapsed_process = time.perf_counter() - start_process_sec

        end_wall_sec = time.time()
        elapsed_wall = end_wall_sec - start_wall_sec
        print("\n{} > Elapsed time {:1.2f} sec ({}); {:1.2f} sec in process".format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_wall_sec)),
            elapsed_wall,
            datetime.timedelta(seconds=elapsed_wall),
            elapsed_process,
        ))


class TestScript(ScriptBase):
    """To test out the command line parser."""

    def define_parameters(self):
        super(TestScript, self).define_parameters()
        self.parser.add_argument('--seed', default='000001', help='simulation seed')
        self.define_option(name="thread", default=1, help="hahahahe", short_name="t")

    def run(self, args):
        print("[TEST] Run args:", args)


if __name__ == '__main__':
    script = TestScript()
    script.cli()
