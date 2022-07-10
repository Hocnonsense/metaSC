# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-06 18:01:40
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-10 12:52:24
 * @FilePath: /metaSC/PyLib/PyLibTool/demo_exec.py
 * @Description:
"""

from datetime import datetime
import logging
import argparse
from typing import Tuple


logger = logging.getLogger(__name__)


def main():
    return 0


def get_args() -> Tuple:
    parser = argparse.ArgumentParser(description=__doc__)
    set_args(parser)
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel.upper())  # info

    parser.print_help()

    return ()


def set_args(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--loglevel", default="INFO", type=str, help="set level of logger"
    )


def run():
    args = get_args()

    now = datetime.now()
    logger.warning(">>> job start at " + now.strftime("%Y-%m-%d %H:%M:%S"))
    state = main(*args)
    logger.warning(">>> job run time: " + str(datetime.now() - now))
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
