# -*- coding: utf-8 -*-
"""
 * @Date: 2020-08-20 16:05:42
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 14:23:26
 * @FilePath: /HScripts/Python/mylib/tool/dict.py
 * @Description: 各种数据结构
        Miscellaneous data structure utilities.
"""


from typing import Iterable, Mapping


def dissoc(mapping, keys):
    # type: (Mapping, Iterable) -> dict
    """Dissociate: Return a new dict like `mapping` without the given `keys`.
    See also `dissoc_strict()`.
    """
    result = dict(mapping)

    for key in keys:
        result.pop(key, None)
    return result


def dissoc_strict(mapping, keys):
    # type: (Mapping, Iterable) -> dict
    """Dissociate: Return a new dict like `mapping` without the given `keys`.
    Raises a KeyError if any of these keys are not in the given mapping.
    """
    result = dict(mapping)

    for key in keys:
        del result[key]
    return result


def select_keys(mapping, keys):
    # type: (Mapping, Iterable) -> dict
    """Return a dict of the entries in mapping with the given keys."""
    return {key: mapping[key] for key in keys}


if __name__ == "__main__":
    dic = {1: 2, 3: 4, 5: 6}
    print("dic:", dic)
    dic1 = select_keys(dic, [1, 3])
    print("select_keys(dic, [2, 4]):", dic1)
    dic2 = dissoc(dic, [1, 2])
    print("dissoc(dic, [1, 2]):", dic2)
    try:
        dic3 = dissoc_strict(dic, [1, 2])
    except KeyError as e:
        print("dissoc_strict(dic, [1, 2])", "raise KeyError:", e)
