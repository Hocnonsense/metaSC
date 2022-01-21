# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-06 18:01:40
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-21 15:37:34
 * @FilePath: /metaSC/PyLib/PyLibTool/demo_exec.1.py
 * @Description:
"""

from datetime import datetime
import click
from typing import Tuple


try:
    from PyLib.PyLibTool.file_info import verbose_import
    logger = verbose_import(__name__, __doc__)
except ImportError:
    import logging
    logger = logging.getLogger(__name__)


def main():
    return 0


@click.command()
@click.option('--loglevel', default='INFO', type=str,
              help='set level of logger')
def get_args(loglevel: str) -> Tuple:
    print(loglevel)
    logging.basicConfig(level=loglevel.upper())  # info

    return ()


def run():
    args = get_args()

    now = datetime.now()
    logger.warning('>>> job start at ' + now.strftime('%Y-%m-%d %H:%M:%S'))
    state = main(*args)
    logger.warning('>>> job run time: ' + str(datetime.now() - now))
    if state == 0:
        logger.info('success!')


if __name__ == '__main__':
    run()
