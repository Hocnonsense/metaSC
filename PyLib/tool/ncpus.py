# -*- coding: utf-8 -*-
"""
 * @Date: 2021-05-18 17:53:35
 * @Editor: Luis Pedro Coelho <luispedro@big-data-biology.org>
 * @LastEditors: Luis Pedro Coelho
 * @LastEditTime: 2022-03-03 16:55:57
 * @FilePath: /metaSC/PyLib/tool/ncpus.py
 * @Description:
"""


def get_ncpus():
    from os import environ

    for ev in ["OMP_NUM_THREADS", "Q_CORES", "Q_CORE"]:
        if ev in environ:
            return int(environ[ev].strip())
    for ev in ["LSB_MCPU_HOSTS"]:
        if ev in environ:
            break
    else:
        return 1
    tokens = environ[ev].strip().split()
    if len(tokens) > 2:
        raise SystemError(
            "Cannot handle this type of environment ({}='{}')".format(ev, environ[ev])
        )
    return int(tokens[1])


if __name__ == "__main__":
    print("Running with {} CPUS.".format(get_ncpus()))
