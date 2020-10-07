# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-07 18:32:46
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-07 18:33:06
 * @FilePath: /HScripts/Python/mylib/tool/table.py
 * @Description:
"""

def pprint_list(list_2d: list, sep: str = " | ") -> str:
    """
     * @description:        print a table with readable seperators
     * @param {list_2d}     a 2-diamond table
                            ! Notice: it don't check the table
     * @return              a readble str
    """
    nrows, ncols = len(list_2d), len(list_2d[0])
    rows = [[] for _ in range(nrows)]
    for c in range(ncols):
        alignment = "^"
        if type(list_2d[1][c]) in [int, float]:
            alignment = ">"
        elif type(list_2d[1][c]) in [str]:
            alignment = "<"
        col_values = [str(list_2d[r][c]) for r in range(nrows)]
        length = max(len(c_value) for c_value in col_values)
        for v, new_v in zip(rows, col_values):
            v.append(("{: " + alignment + str(length) + "}").format(new_v))
    return "\n".join([sep.join(v) for v in rows])
